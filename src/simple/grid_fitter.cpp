#include "simple/grid_fitter.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <queue>
#include <limits>

namespace simple {

namespace {

// 拟合单个格子（指定方向 dir）。
GridCell fitOneCell(const SurfaceWrapper& surf, int r, int c,
                    double u0, double u1, double v0, double v1,
                    const GridConfig& cfg, bool planar, ParamDir dir)
{
    GridCell cell;
    cell.row = r;
    cell.col = c;
    cell.u0 = u0; cell.u1 = u1; cell.v0 = v0; cell.v1 = v1;
    if (planar) {
        cell.plane = fitCellPlane(surf, u0, u1, v0, v1, cfg.nUSamples, cfg.nVSamples);
        cell.fitDir = ParamDir::U;
        cell.maxError = cell.plane.maxError;
        cell.rmsError = cell.plane.rmsError;
    } else {
        cell.ruled = fitCellRuled(surf, u0, u1, v0, v1, dir,
                                  cfg.nUSamples, cfg.nVSamples,
                                  cfg.nRibs, cfg.lambda);
        cell.fitDir = dir;
        cell.maxError = cell.ruled.maxError;
        cell.rmsError = cell.ruled.rmsError;
    }
    return cell;
}

// 全局方向试算：对当前网格所有格子分别按 U、V 方向拟合，取误差总和更小者。
ParamDir determineDir(const SurfaceWrapper& surf,
                      const std::vector<double>& uEdges,
                      const std::vector<double>& vEdges,
                      int nRows, int nCols,
                      const GridConfig& cfg)
{
    double sumU = 0.0, sumV = 0.0;
    for (int r = 0; r < nRows; ++r) {
        for (int c = 0; c < nCols; ++c) {
            double u0 = uEdges[c], u1 = uEdges[c + 1];
            double v0 = vEdges[r], v1 = vEdges[r + 1];
            RuledCellFit fu = fitCellRuled(surf, u0, u1, v0, v1, ParamDir::U,
                                           cfg.nUSamples, cfg.nVSamples,
                                           cfg.nRibs, cfg.lambda);
            RuledCellFit fv = fitCellRuled(surf, u0, u1, v0, v1, ParamDir::V,
                                           cfg.nUSamples, cfg.nVSamples,
                                           cfg.nRibs, cfg.lambda);
            sumU += fu.maxError;
            sumV += fv.maxError;
        }
    }
    return (sumV <= sumU) ? ParamDir::V : ParamDir::U;
}

// 用指定方向重新拟合网格所有格子。
void fitAllCells(GridResult& g, const SurfaceWrapper& surf,
                 const GridConfig& cfg, bool planar, ParamDir dir)
{
    for (int r = 0; r < g.nRows; ++r)
        for (int c = 0; c < g.nCols; ++c)
            g.cells[r * g.nCols + c] = fitOneCell(surf, r, c,
                g.uEdges[c], g.uEdges[c + 1], g.vEdges[r], g.vEdges[r + 1],
                cfg, planar, dir);
}

// 整列切分（在第 c 列插入竖直分割线）：拟合该列所有行的两个子格。
double buildColumnSplit(const SurfaceWrapper& surf, const GridResult& g, int c,
                        const GridConfig& cfg, bool planar, ParamDir dir,
                        std::vector<GridCell>& out)
{
    double um = (g.uEdges[c] + g.uEdges[c + 1]) * 0.5;
    out.clear();
    out.reserve(2 * g.nRows);
    double metric = 0.0;
    for (int r = 0; r < g.nRows; ++r) {
        double v0 = g.vEdges[r], v1 = g.vEdges[r + 1];
        GridCell a = fitOneCell(surf, r, c,     g.uEdges[c], um, v0, v1, cfg, planar, dir);
        GridCell b = fitOneCell(surf, r, c + 1, um, g.uEdges[c + 1], v0, v1, cfg, planar, dir);
        metric = std::max(metric, std::max(a.maxError, b.maxError));
        out.push_back(a);
        out.push_back(b);
    }
    return metric;
}

// 整行切分（在第 r 行插入水平分割线）：拟合该行所有列的两个子格。
double buildRowSplit(const SurfaceWrapper& surf, const GridResult& g, int r,
                     const GridConfig& cfg, bool planar, ParamDir dir,
                     std::vector<GridCell>& out)
{
    double vm = (g.vEdges[r] + g.vEdges[r + 1]) * 0.5;
    out.clear();
    out.reserve(2 * g.nCols);
    double metric = 0.0;
    for (int c = 0; c < g.nCols; ++c) {
        double u0 = g.uEdges[c], u1 = g.uEdges[c + 1];
        GridCell a = fitOneCell(surf, r,     c, u0, u1, g.vEdges[r], vm, cfg, planar, dir);
        GridCell b = fitOneCell(surf, r + 1, c, u0, u1, vm, g.vEdges[r + 1], cfg, planar, dir);
        metric = std::max(metric, std::max(a.maxError, b.maxError));
        out.push_back(a);
        out.push_back(b);
    }
    return metric;
}

void applyColumnSplit(GridResult& g, int c, const std::vector<GridCell>& newCells)
{
    std::vector<double> nu = g.uEdges;
    nu.insert(nu.begin() + c + 1, (g.uEdges[c] + g.uEdges[c + 1]) * 0.5);
    int newCols = g.nCols + 1;
    std::vector<GridCell> nc(g.nRows * newCols);
    for (int r = 0; r < g.nRows; ++r) {
        for (int cc = 0; cc < g.nCols; ++cc) {
            if (cc < c) {
                nc[r * newCols + cc] = g.cells[r * g.nCols + cc];
            } else if (cc == c) {
                nc[r * newCols + c]     = newCells[r * 2 + 0];
                nc[r * newCols + c + 1] = newCells[r * 2 + 1];
            } else {
                nc[r * newCols + cc + 1] = g.cells[r * g.nCols + cc];
                nc[r * newCols + cc + 1].col = cc + 1;
            }
        }
    }
    g.uEdges = std::move(nu);
    g.nCols = newCols;
    g.cells = std::move(nc);
}

void applyRowSplit(GridResult& g, int r, const std::vector<GridCell>& newCells)
{
    std::vector<double> nv = g.vEdges;
    nv.insert(nv.begin() + r + 1, (g.vEdges[r] + g.vEdges[r + 1]) * 0.5);
    int newRows = g.nRows + 1;
    std::vector<GridCell> nc(newRows * g.nCols);
    for (int rr = 0; rr < g.nRows; ++rr) {
        for (int c = 0; c < g.nCols; ++c) {
            if (rr < r) {
                nc[rr * g.nCols + c] = g.cells[rr * g.nCols + c];
            } else if (rr == r) {
                nc[r * g.nCols + c]         = newCells[c * 2 + 0];
                nc[(r + 1) * g.nCols + c]   = newCells[c * 2 + 1];
            } else {
                nc[(rr + 1) * g.nCols + c] = g.cells[rr * g.nCols + c];
                nc[(rr + 1) * g.nCols + c].row = rr + 1;
            }
        }
    }
    g.vEdges = std::move(nv);
    g.nRows = newRows;
    g.cells = std::move(nc);
}

GridResult fitGridImpl(const SurfaceWrapper& surf, const GridConfig& cfg,
                       const std::string& name, bool planar)
{
    GridResult g;
    g.name = name;

    auto [u0, u1] = surf.paramDomainU();
    auto [v0, v1] = surf.paramDomainV();

    int nCols = std::max(1, cfg.nSplitU);
    int nRows = std::max(1, cfg.nSplitV);
    g.nCols = nCols;
    g.nRows = nRows;
    g.uEdges.resize(nCols + 1);
    g.vEdges.resize(nRows + 1);
    for (int c = 0; c <= nCols; ++c) g.uEdges[c] = u0 + (u1 - u0) * c / nCols;
    for (int r = 0; r <= nRows; ++r) g.vEdges[r] = v0 + (v1 - v0) * r / nRows;

    // 初始方向：直纹面由全局试算决定，平面恒为 U
    ParamDir gDir = planar ? ParamDir::U
                           : determineDir(surf, g.uEdges, g.vEdges, nRows, nCols, cfg);
    g.fitDir = gDir;
    g.cells.resize(nRows * nCols);
    fitAllCells(g, surf, cfg, planar, gDir);

    double minULen = (u1 - u0) * 1e-3;
    double minVLen = (v1 - v0) * 1e-3;

    struct Entry { double err; int r, c; bool operator<(const Entry& o) const { return err < o.err; } };

    auto rebuild = [&](std::priority_queue<Entry>& pq) {
        pq = std::priority_queue<Entry>();
        for (int r = 0; r < g.nRows; ++r)
            for (int c = 0; c < g.nCols; ++c) {
                double e = g.cells[r * g.nCols + c].maxError;
                if (e > cfg.tolerance) pq.push({e, r, c});
            }
    };

    std::priority_queue<Entry> pq;
    rebuild(pq);

    int steps = 0;
    while (!pq.empty() && steps < cfg.maxDepth && (int)g.cells.size() < cfg.maxCells) {
        ++steps;
        Entry top = pq.top();
        pq.pop();
        int r = top.r, c = top.c;

        GridCell& cur = g.cells[r * g.nCols + c];
        if (cur.maxError <= cfg.tolerance) continue;

        double du = g.uEdges[c + 1] - g.uEdges[c];
        double dv = g.vEdges[r + 1] - g.vEdges[r];
        bool canU = du > minULen;
        bool canV = dv > minVLen;
        if (!canU && !canV) continue;

        std::vector<GridCell> candU, candV;
        double mU = std::numeric_limits<double>::infinity();
        double mV = std::numeric_limits<double>::infinity();
        if (canU) mU = buildColumnSplit(surf, g, c, cfg, planar, gDir, candU);
        if (canV) mV = buildRowSplit(surf, g, r, cfg, planar, gDir, candV);

        if (mU <= mV) applyColumnSplit(g, c, candU);
        else          applyRowSplit(g, r, candV);

        // 每次二分后重新试算方向（非固定值）
        if (!planar) {
            ParamDir nd = determineDir(surf, g.uEdges, g.vEdges, g.nRows, g.nCols, cfg);
            if (nd != gDir) {
                gDir = nd;
                g.fitDir = nd;
                fitAllCells(g, surf, cfg, planar, gDir);
            }
        }

        rebuild(pq);
    }

    g.maxError = 0.0;
    for (const auto& cell : g.cells)
        g.maxError = std::max(g.maxError, cell.maxError);
    g.toleranceMet = (g.maxError <= cfg.tolerance);

    return g;
}

std::string jsonSafeStr(std::string s) {
    for (char& ch : s) { if (ch == '\\') ch = '/'; if (ch == '"') ch = '\''; }
    return s;
}

} // anon

GridResult fitGridRuled(const SurfaceWrapper& surf, const GridConfig& cfg,
                        const std::string& name)
{
    return fitGridImpl(surf, cfg, name, false);
}

GridResult fitGridPlanar(const SurfaceWrapper& surf, const GridConfig& cfg,
                         const std::string& name)
{
    return fitGridImpl(surf, cfg, name, true);
}

bool exportGridOBJs(const std::string& outDir, const std::string& prefix,
                    const GridResult& gr, bool planar)
{
    for (const auto& cell : gr.cells) {
        int idx = cell.row * gr.nCols + cell.col;
        if (planar) {
            std::string op = outDir + "/" + prefix + "_cell" + std::to_string(idx) + ".obj";
            exportOBJ(op, cell.plane.meshVerts, cell.plane.meshFaces);
            std::string dp = outDir + "/" + prefix + "_cell" + std::to_string(idx) + "_desc.txt";
            std::ofstream o(dp);
            if (o) {
                o << std::fixed << std::setprecision(6);
                o << "centroid = " << cell.plane.centroid.x() << " " << cell.plane.centroid.y()
                  << " " << cell.plane.centroid.z() << "\n";
                o << "normal = " << cell.plane.normal.x() << " " << cell.plane.normal.y()
                  << " " << cell.plane.normal.z() << "\n";
            }
        } else {
            std::string op = outDir + "/" + prefix + "_cell" + std::to_string(idx) + ".obj";
            exportOBJ(op, cell.ruled.ruledMeshVerts, cell.ruled.ruledMeshFaces);
            std::string cp = outDir + "/" + prefix + "_cell" + std::to_string(idx) + "_params.txt";
            exportDirectrixTXT(cp, cell.ruled.curveC0Samples, cell.ruled.curveC1Samples);
        }
    }
    return true;
}

bool exportGridLinesVTK(const std::string& path,
                        const SurfaceWrapper& surf,
                        const GridResult& gr,
                        int nSamplesPerLine)
{
    std::ofstream out(path);
    if (!out) return false;
    if (nSamplesPerLine < 2) nSamplesPerLine = 2;

    auto [uMin, uMax] = surf.paramDomainU();
    auto [vMin, vMax] = surf.paramDomainV();

    struct Line { std::vector<Vec3> pts; };
    std::vector<Line> lines;

    for (double vk : gr.vEdges) {
        Line L;
        for (int i = 0; i < nSamplesPerLine; ++i) {
            double u = uMin + (uMax - uMin) * i / (nSamplesPerLine - 1.0);
            L.pts.push_back(surf.evaluate(u, vk));
        }
        lines.push_back(std::move(L));
    }
    for (double uk : gr.uEdges) {
        Line L;
        for (int j = 0; j < nSamplesPerLine; ++j) {
            double v = vMin + (vMax - vMin) * j / (nSamplesPerLine - 1.0);
            L.pts.push_back(surf.evaluate(uk, v));
        }
        lines.push_back(std::move(L));
    }

    int nPts = 0, nCells = 0, cellSize = 0;
    for (const auto& L : lines) {
        nPts += (int)L.pts.size();
        int segs = (int)L.pts.size() - 1;
        nCells += segs;
        cellSize += segs * 3;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "grid split lines\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << nPts << " float\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& L : lines)
        for (const auto& p : L.pts)
            out << p.x() << " " << p.y() << " " << p.z() << "\n";
    out << "LINES " << nCells << " " << cellSize << "\n";
    int base = 0;
    for (const auto& L : lines) {
        for (int i = 0; i < (int)L.pts.size() - 1; ++i)
            out << "2 " << (base + i) << " " << (base + i + 1) << "\n";
        base += (int)L.pts.size();
    }
    return true;
}

std::string buildGridMetaJson(const std::vector<GridResult>& results,
                              const std::string& mode,
                              const std::string& file1,
                              const std::string& file2)
{
    std::ostringstream o;
    o << "{\"mode\":\"" << mode << "\",\"files\":[\"" << jsonSafeStr(file1)
      << "\",\"" << jsonSafeStr(file2) << "\"],\"surfaces\":[";
    for (size_t s = 0; s < results.size(); ++s) {
        if (s) o << ",";
        const auto& gr = results[s];
        o << "{\"name\":\"" << gr.name << "\",\"nRows\":" << gr.nRows
          << ",\"nCols\":" << gr.nCols
          << ",\"fitDir\":\"" << (gr.fitDir == ParamDir::U ? "U" : "V") << "\""
          << ",\"maxError\":" << gr.maxError
          << ",\"toleranceMet\":" << (gr.toleranceMet ? "true" : "false")
          << ",\"cells\":[";
        for (size_t i = 0; i < gr.cells.size(); ++i) {
            if (i) o << ",";
            const auto& c = gr.cells[i];
            o << "{\"index\":" << (c.row * gr.nCols + c.col)
              << ",\"row\":" << c.row << ",\"col\":" << c.col
              << ",\"u0\":" << c.u0 << ",\"u1\":" << c.u1
              << ",\"v0\":" << c.v0 << ",\"v1\":" << c.v1
              << ",\"fitDir\":\"" << (c.fitDir == ParamDir::U ? "U" : "V") << "\""
              << ",\"maxErr\":" << c.maxError
              << ",\"rmsErr\":" << c.rmsError << "}";
        }
        o << "]}";
    }
    o << "]}";
    return o.str();
}

} // namespace simple
