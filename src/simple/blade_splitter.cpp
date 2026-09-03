#include "simple/blade_splitter.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace simple
{

    BladeSplitResult splitBladeFaceBySection(
        const TopoDS_Face &face, int nSec, int nPts, int hw, double, double)
    {
        BladeSplitResult result;
        std::ostringstream log;

        BRepAdaptor_Surface adapt(face, true);
        Handle(Geom_BSplineSurface) surf = adapt.BSpline();
        if (surf.IsNull())
            surf = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(face));
        if (surf.IsNull())
        {
            result.message = "No BSpline";
            return result;
        }

        double uMin = adapt.FirstUParameter(), uMax = adapt.LastUParameter();
        double vMin = adapt.FirstVParameter(), vMax = adapt.LastVParameter();
        double uRng = uMax - uMin, vRng = vMax - vMin;
        if (uRng < 1e-6)
        {
            result.message = "U too small";
            return result;
        }

        nSec = std::max(3, nSec);
        nPts = std::max(50, nPts);
        double um = uRng * 0.08;
        std::vector<double> uH;
        for (int i = 0; i < nSec; ++i)
            uH.push_back(uMin + um + (uRng - 2 * um) * i / (nSec - 1));

        std::vector<double> cvAvg(nPts, 0.0);
        std::vector<double> cvMax(nPts, 0.0);
        std::vector<int> ct(nPts, 0);
        for (double uh : uH)
        {
            Handle(Geom_Curve) c = surf->UIso(uh);
            if (c.IsNull())
                continue;
            GeomAdaptor_Curve gac(c);
            double t0 = gac.FirstParameter(), t1 = gac.LastParameter();
            for (int i = 0; i < nPts; ++i)
            {
                double t = t0 + (t1 - t0) * i / (nPts - 1);
                gp_Pnt p;
                gp_Vec d1, d2;
                gac.D2(t, p, d1, d2);
                double m = d1.Magnitude();
                double k = (m > 1e-10) ? d1.Crossed(d2).Magnitude() / (m * m * m) : 0.0;
                cvAvg[i] += k;
                cvMax[i] = std::max(cvMax[i], k);
                ct[i]++;
            }
        }
        for (int i = 0; i < nPts; ++i)
            if (ct[i] > 0)
                cvAvg[i] /= ct[i];

        // 切分信号：同时考虑最值与平均值（取更大者），保留局部卷曲角的高曲率
        std::vector<double> cv(nPts);
        for (int i = 0; i < nPts; ++i)
            cv[i] = std::max(cvAvg[i], cvMax[i]);

        std::vector<double> sm(nPts);
        int h = 1;
        for (int i = 0; i < nPts; ++i)
        {
            double s = 0;
            int c = 0;
            for (int j = std::max(0, i - h); j < std::min(nPts, i + h + 1); ++j)
            {
                s += cv[j];
                ++c;
            }
            sm[i] = c > 0 ? s / c : 0;
        }

        double mx = *std::max_element(sm.begin(), sm.end());
        double mn = *std::min_element(sm.begin(), sm.end());
        log << " maxC=" << mx;

        // —— 叶缘识别：叶缘 = 曲率谷（窄，远小于叶盆叶背）。
        // 切分点 = 曲率急剧下降的起点（进入谷）与急剧上升的终点（离开谷），而非峰谷中点。
        auto vOf = [&](int i)
        { return vMin + vRng * (double)i / (nPts - 1.0); };

        double thDrop = (mx - mn) * 0.05; // 阈值调小 → 更早识别曲率下降，切除范围更大
        std::vector<double> bp;
        for (int i = 0; i < nPts; ++i)
        {
            double prv = sm[(i - 1 + nPts) % nPts];
            double cur = sm[i];
            double nxt = sm[(i + 1) % nPts];
            if (prv - cur > thDrop) // 骤降：谷起点
                bp.push_back(vOf(i));
            if (nxt - cur > thDrop) // 骤升：谷终点
                bp.push_back(vOf((i + 1) % nPts));
        }
        std::sort(bp.begin(), bp.end());
        std::vector<double> bpu;
        for (double v : bp)
            if (bpu.empty() || v - bpu.back() > vRng * 0.005)
                bpu.push_back(v);
        bp = bpu;

        auto rcWrap = [&](double a, double b)
        {
            int ia = (int)((a - vMin) / vRng * (nPts - 1.0));
            int ib = (int)((b - vMin) / vRng * (nPts - 1.0));
            double s = 0;
            int c = 0;
            for (int i = ia;; ++i)
            {
                s += sm[i % nPts];
                ++c;
                if (i % nPts == ib % nPts)
                    break;
            }
            return c > 0 ? s / c : 0.0;
        };

        // 构建区域（环绕）：相邻切分点之间为一段
        std::vector<SplitRegion> regs;
        int nCuts = (int)bp.size();
        if (nCuts >= 2)
        {
            for (int k = 0; k < nCuts; ++k)
            {
                double a = bp[k];
                double b = bp[(k + 1) % nCuts];
                double end = (k == nCuts - 1) ? b + vRng : b; // 最后一段跨环绕
                SplitRegion r;
                r.uStart = a;
                r.uEnd = end;
                r.avgCurv = rcWrap(a, end);
                regs.push_back(r);
            }
        }
        else
        {
            SplitRegion r;
            r.uStart = vMin;
            r.uEnd = vMax;
            r.avgCurv = rcWrap(vMin, vMax);
            regs.push_back(r);
        }

        double rMax = 0, rMin = 1e30;
        for (auto &r : regs)
        {
            if (r.avgCurv > rMax)
                rMax = r.avgCurv;
            if (r.avgCurv < rMin)
                rMin = r.avgCurv;
        }
        double rMed = (rMax + rMin) * 0.5;

        for (auto &r : regs)
        {
            double L = r.uEnd - r.uStart;
            if (L < vRng * 0.08)
                r.label = "edge";
            else if (r.avgCurv >= rMed)
                r.label = "suction";
            else
                r.label = "pressure";
        }

        std::vector<SplitRegion> mr;
        for (auto &r : regs)
        {
            if (!mr.empty() && mr.back().label == r.label)
            {
                mr.back().uEnd = r.uEnd;
                mr.back().avgCurv = (mr.back().avgCurv + r.avgCurv) * 0.5;
            }
            else
                mr.push_back(r);
        }

        if (mr.size() == 3 && mr[1].label == "pressure" && mr[1].uEnd - mr[1].uStart > vRng * 0.6)
        {
            SplitRegion side = mr[1];
            mr.erase(mr.begin() + 1);
            double mid = (side.uStart + side.uEnd) * 0.5;
            SplitRegion r1, r2;
            r1.uStart = side.uStart;
            r1.uEnd = mid;
            r2.uStart = mid;
            r2.uEnd = side.uEnd;
            int is = std::max(0, (int)((side.uStart - vMin) / vRng * (nPts - 1)));
            int im = std::max(0, (int)((mid - vMin) / vRng * (nPts - 1)));
            int ie = std::min(nPts - 1, (int)((side.uEnd - vMin) / vRng * (nPts - 1)));
            double s1 = 0, s2 = 0;
            int c1 = 0, c2 = 0;
            for (int i = is; i < im; i++)
            {
                s1 += sm[i];
                c1++;
            }
            for (int i = im; i <= ie; i++)
            {
                s2 += sm[i];
                c2++;
            }
            r1.avgCurv = c1 > 0 ? s1 / c1 : 0;
            r2.avgCurv = c2 > 0 ? s2 / c2 : 0;
            r1.label = (r1.avgCurv >= r2.avgCurv) ? "suction" : "pressure";
            r2.label = (r1.label == "suction") ? "pressure" : "suction";
            mr.insert(mr.begin() + 1, r1);
            mr.insert(mr.begin() + 2, r2);
        }

        result.regions = mr;
        result.success = true;
        result.splitDir = 'V';
        for (auto &r : mr)
            log << "\n  V[" << r.uStart << "," << r.uEnd << "] L="
                << (r.uEnd - r.uStart) << " C=" << r.avgCurv << " " << r.label;
        result.message = log.str();
        return result;
    }

} // namespace simple
