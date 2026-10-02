#include "simple/machining.hpp"
#include "simple/common.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <chrono>
#include <vector>
#include <array>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <functional>

namespace simple {
namespace {

namespace fs = std::filesystem;
using Face = std::array<int, 3>;

// ────────────────────────────────────────────────────────────
// 最小 JSON 解析器（仅用于读取 meta.json）
// ────────────────────────────────────────────────────────────
struct Json {
    enum Type { Null, Bool, Number, String, Array, Object };
    Type type = Null;
    bool b = false;
    double num = 0.0;
    std::string str;
    std::vector<Json> arr;
    std::vector<std::pair<std::string, Json>> obj;

    const Json* find(const std::string& key) const {
        for (const auto& kv : obj)
            if (kv.first == key) return &kv.second;
        return nullptr;
    }
    int asInt() const { return (int)std::lround(num); }
    double asDouble() const { return num; }
    const std::string& asStr() const { return str; }
};

class JsonParser {
public:
    explicit JsonParser(const std::string& s) : s_(s), i_(0) {}
    bool parse(Json& out) {
        skipWs();
        return parseValue(out);
    }

private:
    std::string s_;
    size_t i_;

    void skipWs() {
        while (i_ < s_.size() && std::isspace((unsigned char)s_[i_])) ++i_;
    }
    bool parseValue(Json& out) {
        skipWs();
        if (i_ >= s_.size()) return false;
        char c = s_[i_];
        if (c == '{') return parseObject(out);
        if (c == '[') return parseArray(out);
        if (c == '"') return parseString(out);
        if (c == 't' || c == 'f') return parseBool(out);
        if (c == 'n') return parseNull(out);
        return parseNumber(out);
    }
    bool parseObject(Json& out) {
        out.type = Json::Object;
        ++i_;  // {
        skipWs();
        if (i_ < s_.size() && s_[i_] == '}') { ++i_; return true; }
        while (true) {
            skipWs();
            Json key;
            if (!parseString(key)) return false;
            skipWs();
            if (i_ >= s_.size() || s_[i_] != ':') return false;
            ++i_;
            Json val;
            if (!parseValue(val)) return false;
            out.obj.push_back({key.str, std::move(val)});
            skipWs();
            if (i_ < s_.size() && s_[i_] == ',') { ++i_; continue; }
            if (i_ < s_.size() && s_[i_] == '}') { ++i_; return true; }
            return false;
        }
    }
    bool parseArray(Json& out) {
        out.type = Json::Array;
        ++i_;  // [
        skipWs();
        if (i_ < s_.size() && s_[i_] == ']') { ++i_; return true; }
        while (true) {
            Json val;
            if (!parseValue(val)) return false;
            out.arr.push_back(std::move(val));
            skipWs();
            if (i_ < s_.size() && s_[i_] == ',') { ++i_; continue; }
            if (i_ < s_.size() && s_[i_] == ']') { ++i_; return true; }
            return false;
        }
    }
    bool parseString(Json& out) {
        out.type = Json::String;
        ++i_;  // "
        out.str.clear();
        while (i_ < s_.size()) {
            char c = s_[i_];
            if (c == '"') { ++i_; return true; }
            if (c == '\\' && i_ + 1 < s_.size()) {
                ++i_;
                char e = s_[i_];
                if (e == '"') out.str += '"';
                else if (e == '\\') out.str += '\\';
                else if (e == '/') out.str += '/';
                else if (e == 'n') out.str += '\n';
                else if (e == 't') out.str += '\t';
                else if (e == 'r') out.str += '\r';
                else out.str += e;
                ++i_;
            } else {
                out.str += c;
                ++i_;
            }
        }
        return false;
    }
    bool parseBool(Json& out) {
        out.type = Json::Bool;
        if (s_.compare(i_, 4, "true") == 0) { out.b = true; i_ += 4; return true; }
        if (s_.compare(i_, 5, "false") == 0) { out.b = false; i_ += 5; return true; }
        return false;
    }
    bool parseNull(Json& out) {
        out.type = Json::Null;
        if (s_.compare(i_, 4, "null") == 0) { i_ += 4; return true; }
        return false;
    }
    bool parseNumber(Json& out) {
        out.type = Json::Number;
        size_t start = i_;
        if (i_ < s_.size() && (s_[i_] == '-' || s_[i_] == '+')) ++i_;
        while (i_ < s_.size() &&
               (std::isdigit((unsigned char)s_[i_]) || s_[i_] == '.' ||
                s_[i_] == 'e' || s_[i_] == 'E' || s_[i_] == '-' || s_[i_] == '+'))
            ++i_;
        if (i_ == start) return false;
        try {
            out.num = std::stod(s_.substr(start, i_ - start));
        } catch (...) {
            return false;
        }
        return true;
    }
};

// ────────────────────────────────────────────────────────────
// 向量工具（Eigen Vec3）
// ────────────────────────────────────────────────────────────
inline Vec3 normalized(const Vec3& a) {
    double m = a.norm();
    return m < 1e-12 ? Vec3(0, 0, 0) : a / m;
}
inline Vec3 lerp3(const Vec3& a, const Vec3& b, double t) { return a + (b - a) * t; }
inline double angleDeg(const Vec3& a, const Vec3& b) {
    double m = a.norm() * b.norm();
    if (m < 1e-12) return 0.0;
    double d = clamp(a.dot(b) / m, -1.0, 1.0);
    return radToDeg(std::acos(d));
}
inline Vec3Arr tangents(const Vec3Arr& C) {
    int n = (int)C.size();
    Vec3Arr out(n, Vec3(0, 0, 0));
    for (int i = 0; i < n; ++i) {
        if (n < 2) out[i] = Vec3(0, 0, 0);
        else if (i == 0) out[i] = normalized(C[1] - C[0]);
        else if (i == n - 1) out[i] = normalized(C[n - 1] - C[n - 2]);
        else out[i] = normalized(C[i + 1] - C[i - 1]);
    }
    return out;
}

// ────────────────────────────────────────────────────────────
// 文件读取
// ────────────────────────────────────────────────────────────
std::pair<Vec3Arr, Vec3Arr> loadDirectrix(const std::string& path) {
    Vec3Arr C0, C1;
    std::ifstream f(path);
    if (!f) return {C0, C1};
    std::string line, cur;
    while (std::getline(f, line)) {
        // 去掉首尾空白
        size_t a = line.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) continue;
        size_t b = line.find_last_not_of(" \t\r\n");
        line = line.substr(a, b - a + 1);
        if (line.empty()) continue;
        if (line[0] == '[') {
            cur = line.substr(1, line.find(']') - 1);
            continue;
        }
        if ((cur == "C0" || cur == "C1") &&
            line[0] != 'n' && line[0] != 'd' && line[0] != 'm') {
            std::istringstream ss(line);
            double x, y, z;
            if (ss >> x >> y >> z) {
                (cur == "C0" ? C0 : C1).push_back(Vec3(x, y, z));
            }
        }
    }
    return {C0, C1};
}

std::pair<Vec3Arr, std::vector<Face>> loadObj(const std::string& path) {
    Vec3Arr verts;
    std::vector<Face> faces;
    std::ifstream f(path);
    if (!f) return {verts, faces};
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        std::string tok;
        ss >> tok;
        if (tok == "v") {
            double x, y, z;
            if (ss >> x >> y >> z) verts.push_back(Vec3(x, y, z));
        } else if (tok == "f") {
            std::vector<int> idxs;
            std::string p;
            while (ss >> p) {
                size_t slash = p.find('/');
                std::string si = (slash == std::string::npos) ? p : p.substr(0, slash);
                if (!si.empty()) {
                    try { idxs.push_back(std::stoi(si) - 1); } catch (...) {}
                }
            }
            if (idxs.size() == 3) faces.push_back({idxs[0], idxs[1], idxs[2]});
        }
    }
    return {verts, faces};
}

double meshArea(const Vec3Arr& verts, const std::vector<Face>& faces) {
    double area = 0.0;
    for (const auto& f : faces) {
        const Vec3& a = verts[f[0]];
        const Vec3& b = verts[f[1]];
        const Vec3& c = verts[f[2]];
        area += 0.5 * (b - a).cross(c - a).norm();
    }
    return area;
}

double readOriginalArea(const std::string& inputDir) {
    double total = 0.0;
    bool found = false;
    for (const auto& e : fs::directory_iterator(inputDir)) {
        if (!e.is_regular_file()) continue;
        std::string name = e.path().filename().string();
        if (name.rfind("blade", 0) == 0 && name.size() >= 9 &&
            name.substr(name.size() - 9) == "_mesh.obj") {
            auto [verts, faces] = loadObj(e.path().string());
            if (!verts.empty() && !faces.empty()) {
                total += meshArea(verts, faces);
                found = true;
            }
        }
    }
    return found ? total : -1.0;  // -1 表示未找到
}

// ────────────────────────────────────────────────────────────
// 几何量
// ────────────────────────────────────────────────────────────
double curveLength(const Vec3Arr& C) {
    double L = 0.0;
    for (size_t i = 0; i + 1 < C.size(); ++i)
        L += (C[i + 1] - C[i]).norm();
    return L;
}

double meanRulingLength(const Vec3Arr& C0, const Vec3Arr& C1) {
    int n = (int)std::min(C0.size(), C1.size());
    if (n == 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += (C1[i] - C0[i]).norm();
    return s / n;
}

struct RulingNormals {
    Vec3Arr n0, n1, r;
};

RulingNormals rulingNormals(const Vec3Arr& C0, const Vec3Arr& C1) {
    int n = (int)std::min(C0.size(), C1.size());
    RulingNormals out;
    out.n0.resize(n, Vec3(0, 0, 0));
    out.n1.resize(n, Vec3(0, 0, 0));
    out.r.resize(n, Vec3(0, 0, 0));
    Vec3Arr T0 = tangents(C0);
    Vec3Arr T1 = tangents(C1);
    for (int i = 0; i < n; ++i) {
        Vec3 r = C1[i] - C0[i];
        out.r[i] = r;
        out.n0[i] = normalized(T0[i].cross(r));
        out.n1[i] = normalized(T1[i].cross(r));
    }
    return out;
}

double twistAngleDeg(const Vec3Arr& C0, const Vec3Arr& C1) {
    RulingNormals rn = rulingNormals(C0, C1);
    double mx = 0.0;
    for (size_t i = 0; i < rn.n0.size(); ++i)
        mx = std::max(mx, angleDeg(rn.n0[i], rn.n1[i]));
    return mx;
}

double surfaceArea(const Vec3Arr& C0, const Vec3Arr& C1) {
    int n = (int)std::min(C0.size(), C1.size());
    double area = 0.0;
    for (int i = 0; i + 1 < n; ++i) {
        Vec3 r = C1[i] - C0[i];
        Vec3 dC = C0[i + 1] - C0[i];
        area += r.cross(dC).norm();
    }
    return area;
}

// ────────────────────────────────────────────────────────────
// 加工时间模型
// ────────────────────────────────────────────────────────────
double flankCutTime(const Vec3Arr& C0, double feed) {
    return curveLength(C0) / feed * 60.0;
}

double pointStepover(double ballR, double scallop) {
    return 2.0 * std::sqrt(std::max(0.0, 2.0 * ballR * scallop - scallop * scallop));
}

double pointCutTime(const Vec3Arr& C0, const Vec3Arr& C1, double feed,
                    double ballR, double scallop) {
    double step = pointStepover(ballR, scallop);
    double A = surfaceArea(C0, C1);
    double totalLen = step > 0 ? A / step : 0.0;
    return totalLen / feed * 60.0;
}

// ────────────────────────────────────────────────────────────
// 刀轨（CL）折线
// ────────────────────────────────────────────────────────────
std::vector<int> sampleIndices(int n, int nAlong) {
    std::vector<int> idx;
    if (n <= nAlong) {
        for (int i = 0; i < n; ++i) idx.push_back(i);
    } else {
        for (int i = 0; i < nAlong; ++i)
            idx.push_back((int)std::lround(i * (n - 1) / (double)(nAlong - 1)));
    }
    return idx;
}

// 返回 (feed_pts, axis_segs)
std::pair<Vec3Arr, std::vector<std::array<Vec3, 2>>> flankClLines(
    const Vec3Arr& C0, const Vec3Arr& C1, double toolR, double flip, int nAlong) {
    Vec3Arr feed;
    std::vector<std::array<Vec3, 2>> axes;
    int n = (int)C0.size();
    if (n < 2) return {feed, axes};
    auto idx = sampleIndices(n, nAlong);
    RulingNormals rn = rulingNormals(C0, C1);
    double ext = toolR;
    for (int i : idx) {
        Vec3 mid = lerp3(C0[i], C1[i], 0.5);
        Vec3 nm = normalized(rn.n0[i] + rn.n1[i]) * flip;
        Vec3 r = rn.r[i];
        double rlen = r.norm();
        if (rlen < 1e-12) continue;
        Vec3 rhat = r / rlen;
        Vec3 axisP = mid + nm * toolR;
        feed.push_back(axisP);
        double half = rlen * 0.5 + ext;
        axes.push_back({axisP - rhat * half, axisP + rhat * half});
    }
    return {feed, axes};
}

// 返回多条行切折线
std::vector<Vec3Arr> pointClLines(const Vec3Arr& C0, const Vec3Arr& C1,
                                  double stepover, double ballR, double flip,
                                  int nAlong) {
    std::vector<Vec3Arr> lines;
    int n = (int)C0.size();
    if (n < 2) return lines;
    auto idx = sampleIndices(n, nAlong);
    double meanR = meanRulingLength(C0, C1);
    int nAcross = stepover > 0 ? std::max(1, (int)std::lround(meanR / stepover)) : 1;
    RulingNormals rn = rulingNormals(C0, C1);
    for (int j = 0; j <= nAcross; ++j) {
        double t = nAcross > 0 ? (double)j / nAcross : 0.0;
        Vec3Arr row;
        row.reserve(idx.size());
        for (int i : idx) {
            Vec3 nm = normalized(lerp3(rn.n0[i], rn.n1[i], t)) * flip;
            Vec3 surf = lerp3(C0[i], C1[i], t);
            row.push_back(surf + nm * ballR);
        }
        lines.push_back(std::move(row));
    }
    return lines;
}

// ────────────────────────────────────────────────────────────
// 折线导出（VTK / CSV / DXF）
// ────────────────────────────────────────────────────────────
void writeVtkPolylines(const std::string& path,
                       const std::vector<Vec3Arr>& lines) {
    std::ofstream out(path);
    if (!out) return;
    size_t nPts = 0, nSegs = 0;
    for (const auto& l : lines) {
        nPts += l.size();
        if (l.size() >= 2) nSegs += l.size() - 1;
    }
    out << "# vtk DataFile Version 3.0\n";
    out << "toolpath\nASCII\nDATASET POLYDATA\n";
    out << "POINTS " << nPts << " float\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& l : lines)
        for (const auto& p : l)
            out << p.x() << " " << p.y() << " " << p.z() << "\n";
    out << "LINES " << nSegs << " " << nSegs * 3 << "\n";
    int base = 0;
    for (const auto& l : lines) {
        for (size_t k = 0; k + 1 < l.size(); ++k)
            out << "2 " << (base + k) << " " << (base + k + 1) << "\n";
        base += (int)l.size();
    }
}

void writeLinesCsv(const std::string& path, const std::vector<Vec3Arr>& lines) {
    std::ofstream out(path);
    if (!out) return;
    out << "line_id,point_idx,x,y,z\n";
    out << std::fixed << std::setprecision(6);
    for (size_t li = 0; li < lines.size(); ++li)
        for (size_t pi = 0; pi < lines[li].size(); ++pi)
            out << li << "," << pi << ","
                << lines[li][pi].x() << "," << lines[li][pi].y() << ","
                << lines[li][pi].z() << "\n";
}

void writeLinesDxf(const std::string& path, const std::vector<Vec3Arr>& lines) {
    std::ofstream out(path);
    if (!out) return;
    out << "0\nSECTION\n2\nENTITIES\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& l : lines) {
        if (l.size() < 2) continue;
        out << "0\nPOLYLINE\n8\n0\n66\n1\n70\n8\n";
        for (const auto& p : l) {
            out << "0\nVERTEX\n8\n0\n";
            out << "10\n" << p.x() << "\n20\n" << p.y() << "\n30\n" << p.z() << "\n";
        }
        out << "0\nSEQEND\n";
    }
    out << "0\nENDSEC\n0\nEOF\n";
}

// ────────────────────────────────────────────────────────────
// 数据模型与连通域
// ────────────────────────────────────────────────────────────
struct Patch {
    std::string name;
    int blade = 0;
    int cellIdx = -1;
    int row = -1, col = -1;
    Vec3Arr C0, C1;
    double directrixLen = 0.0;
    double meanRuling = 0.0;
    double area = 0.0;
    double twist = 0.0;
    double flankTime = 0.0;
    double pointTime = 0.0;
    double stepover = 0.0;
    double normalFlip = 1.0;
};

std::pair<int, int> parseCellName(const std::string& name) {
    int blade = (name.find("blade1") != std::string::npos) ? 0 : 1;
    int idx = -1;
    size_t pos = name.find("_cell");
    if (pos != std::string::npos) {
        std::string tail = name.substr(pos + 5);
        try { idx = std::stoi(tail); } catch (...) { idx = -1; }
    }
    return {blade, idx};
}

// 4 邻接连通域（并查集），返回连通域数量
int countConnectedRegions(const std::vector<Patch>& patches) {
    if (patches.empty()) return 0;
    // 无网格元数据时退化为每片一个区域
    if (std::any_of(patches.begin(), patches.end(),
                    [](const Patch& p) { return p.row < 0 || p.col < 0; }))
        return (int)patches.size();

    auto key = [](int blade, int row, int col) { return (blade * 1000000) + row * 10000 + col; };
    std::unordered_map<int, int> parent;
    for (const auto& p : patches) parent[key(p.blade, p.row, p.col)] = key(p.blade, p.row, p.col);

    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[ra] = rb;
    };

    for (const auto& p : patches) {
        int k = key(p.blade, p.row, p.col);
        for (auto [dr, dc] : {std::pair{0, 1}, std::pair{1, 0}}) {
            int nk = key(p.blade, p.row + dr, p.col + dc);
            if (parent.count(nk)) unite(k, nk);
        }
    }

    std::unordered_map<int, int> groups;
    for (const auto& p : patches) groups[find(key(p.blade, p.row, p.col))]++;
    return (int)groups.size();
}

// ────────────────────────────────────────────────────────────
// 法向翻转（刀具偏移到曲面外侧）
// ────────────────────────────────────────────────────────────
std::unordered_map<int, double> bladeNormalFlips(const std::vector<Patch>& patches) {
    std::unordered_map<int, Vec3> cents;
    std::unordered_map<int, int> counts;
    for (const auto& p : patches) {
        int b = p.blade;
        if (!cents.count(b)) cents[b] = Vec3(0, 0, 0);
        Vec3& c = cents[b];
        for (const auto& v : p.C0) c += v;
        for (const auto& v : p.C1) c += v;
        counts[b] += (int)(p.C0.size() + p.C1.size());
    }
    for (auto& [b, c] : cents) {
        if (counts[b] > 0) c /= (double)counts[b];
    }

    std::unordered_map<int, double> flips;
    for (int blade = 0; blade <= 1; ++blade) {
        int other = 1 - blade;
        if (!cents.count(blade) || !cents.count(other)) {
            flips[blade] = 1.0;
            continue;
        }
        Vec3 ref = cents[other];
        double s = 0.0;
        int n = 0;
        for (const auto& p : patches) {
            if (p.blade != blade) continue;
            RulingNormals rn = rulingNormals(p.C0, p.C1);
            for (size_t i = 0; i < rn.n0.size(); ++i) {
                Vec3 mid = lerp3(p.C0[i], p.C1[i], 0.5);
                s += rn.n0[i].dot(mid - ref);
                ++n;
            }
        }
        flips[blade] = (s >= 0.0) ? 1.0 : -1.0;
    }
    return flips;
}

// ────────────────────────────────────────────────────────────
// summary.json
// ────────────────────────────────────────────────────────────
std::string jsonSafe(const std::string& s) {
    std::string r = s;
    for (char& c : r) {
        if (c == '\\') c = '/';
        if (c == '"') c = '\'';
    }
    return r;
}

} // anon

// ============================================================
// 主入口
// ============================================================
MachiningSummary computeToolpath(const std::string& inputDir,
                                 const std::string& outputDir,
                                 const MachiningConfig& cfg) {
    MachiningSummary sum;
    sum.ok = false;

    auto t0 = std::chrono::steady_clock::now();

    fs::create_directories(outputDir);

    // 1) 收集 *_params.txt
    std::vector<std::string> paramFiles;
    for (const auto& e : fs::directory_iterator(inputDir)) {
        if (!e.is_regular_file()) continue;
        std::string name = e.path().filename().string();
        if (name.size() >= 11 && name.substr(name.size() - 11) == "_params.txt")
            paramFiles.push_back(e.path().string());
    }
    if (paramFiles.empty()) {
        sum.errorMsg = "No *_params.txt found in: " + inputDir;
        return sum;
    }
    std::sort(paramFiles.begin(), paramFiles.end());

    // 2) 读 meta.json 得到每面的 (nRows, nCols)
    std::unordered_map<int, std::pair<int, int>> gridMeta;  // blade -> (nRows, nCols)
    {
        std::ifstream f(inputDir + "/meta.json");
        if (f) {
            std::stringstream ss;
            ss << f.rdbuf();
            Json root;
            JsonParser parser(ss.str());
            if (parser.parse(root)) {
                const Json* surfaces = root.find("surfaces");
                if (surfaces && surfaces->type == Json::Array) {
                    for (size_t i = 0; i < surfaces->arr.size(); ++i) {
                        const Json& s = surfaces->arr[i];
                        const Json* nr = s.find("nRows");
                        const Json* nc = s.find("nCols");
                        if (nr && nc)
                            gridMeta[(int)i] = {nr->asInt(), nc->asInt()};
                    }
                }
            }
        }
    }

    // 3) 逐片计算
    std::vector<Patch> patches;
    for (const auto& fp : paramFiles) {
        auto [C0, C1] = loadDirectrix(fp);
        if (C0.size() < 2 || C1.size() < 2) continue;
        Patch p;
        std::string stem = fs::path(fp).filename().string();
        stem = stem.substr(0, stem.size() - 11);  // 去掉 _params.txt
        p.name = stem;
        std::tie(p.blade, p.cellIdx) = parseCellName(stem);
        if (gridMeta.count(p.blade) && gridMeta[p.blade].second > 0 && p.cellIdx >= 0) {
            int nCols = gridMeta[p.blade].second;
            p.row = p.cellIdx / nCols;
            p.col = p.cellIdx % nCols;
        }
        p.C0 = std::move(C0);
        p.C1 = std::move(C1);
        p.directrixLen = curveLength(p.C0);
        p.meanRuling = meanRulingLength(p.C0, p.C1);
        p.area = surfaceArea(p.C0, p.C1);
        p.twist = twistAngleDeg(p.C0, p.C1);
        p.flankTime = flankCutTime(p.C0, cfg.feed);
        p.stepover = pointStepover(cfg.ball_r, cfg.scallop);
        p.pointTime = pointCutTime(p.C0, p.C1, cfg.feed, cfg.ball_r, cfg.scallop);
        patches.push_back(std::move(p));
    }
    if (patches.empty()) {
        sum.errorMsg = "No valid directrices parsed from: " + inputDir;
        return sum;
    }

    // 4) 法向翻转
    auto flips = bladeNormalFlips(patches);
    for (auto& p : patches) p.normalFlip = flips.count(p.blade) ? flips[p.blade] : 1.0;

    // 5) 原始曲面面积（点铣基线）
    double originalArea = readOriginalArea(inputDir);

    // 6) 汇总时间
    int flankRegions = countConnectedRegions(patches);
    double flankCut = 0.0;
    for (const auto& p : patches) flankCut += p.flankTime;
    double flankOverhead = (double)flankRegions * cfg.overhead;
    double flankTotal = flankCut + flankOverhead;

    double stepover = pointStepover(cfg.ball_r, cfg.scallop);
    double pointCut = 0.0;
    if (originalArea > 0.0 && stepover > 0.0)
        pointCut = originalArea / stepover / cfg.feed * 60.0;
    else
        for (const auto& p : patches) pointCut += p.pointTime;
    double pointTotal = pointCut + cfg.point_overhead;

    double totalArea = 0.0;
    for (const auto& p : patches) totalArea += p.area;

    double speedup = (flankTotal > 0.0) ? pointTotal / flankTotal : 0.0;

    // 7) 导出刀轨（侧铣 + 点铣）
    {
        std::vector<Vec3Arr> flankLines, feedLines, axisLines, pointLines;
        for (const auto& p : patches) {
            auto [feed, axes] = flankClLines(p.C0, p.C1, cfg.tool_r, p.normalFlip, 50);
            if (feed.size() >= 2) {
                flankLines.push_back(feed);
                feedLines.push_back(feed);
            }
            for (const auto& seg : axes) {
                flankLines.push_back({seg[0], seg[1]});
                axisLines.push_back({seg[0], seg[1]});
            }
            for (auto& l : pointClLines(p.C0, p.C1, stepover, cfg.ball_r, p.normalFlip, 50))
                pointLines.push_back(std::move(l));
        }
        // 合并（进给 + 刀轴），供 NX/DXF 导入
        writeVtkPolylines(outputDir + "/toolpath_flank.vtk", flankLines);
        writeLinesCsv(outputDir + "/toolpath_flank.csv", flankLines);
        writeLinesDxf(outputDir + "/toolpath_flank.dxf", flankLines);
        // 分开的进给轨迹与刀轴，供 UI 分别可视化
        writeVtkPolylines(outputDir + "/toolpath_flank_feed.vtk", feedLines);
        writeVtkPolylines(outputDir + "/toolpath_flank_axis.vtk", axisLines);
        // 点铣
        writeVtkPolylines(outputDir + "/toolpath_point.vtk", pointLines);
        writeLinesCsv(outputDir + "/toolpath_point.csv", pointLines);
        writeLinesDxf(outputDir + "/toolpath_point.dxf", pointLines);
    }

    // 8) summary.json
    sum.elapsedSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    {
        std::ofstream o(outputDir + "/summary.json");
        if (o) {
            o << std::fixed << std::setprecision(6);
            o << "{\"ok\":true,\"config\":{"
              << "\"feed\":" << cfg.feed
              << ",\"tool_r\":" << cfg.tool_r
              << ",\"ball_r\":" << cfg.ball_r
              << ",\"scallop\":" << cfg.scallop
              << ",\"overhead\":" << cfg.overhead
              << ",\"point_overhead\":" << cfg.point_overhead << "},"
              << "\"num_patches\":" << patches.size()
              << ",\"flank_regions\":" << flankRegions
              << ",\"total_area\":" << totalArea
              << ",\"original_area\":" << (originalArea > 0.0 ? originalArea : 0.0)
              << ",\"flank\":{\"cut\":" << flankCut
              << ",\"overhead\":" << flankOverhead
              << ",\"total\":" << flankTotal << "},"
              << "\"point\":{\"cut\":" << pointCut
              << ",\"overhead\":" << cfg.point_overhead
              << ",\"total\":" << pointTotal << "},"
              << "\"speedup\":" << speedup
              << ",\"elapsed_sec\":" << sum.elapsedSec
              << ",\"patches\":[";
            for (size_t i = 0; i < patches.size(); ++i) {
                if (i) o << ",";
                const auto& p = patches[i];
                o << "{\"name\":\"" << jsonSafe(p.name) << "\""
                  << ",\"blade\":" << p.blade
                  << ",\"twist\":" << p.twist
                  << ",\"directrix_len\":" << p.directrixLen
                  << ",\"area\":" << p.area
                  << ",\"flank_time\":" << p.flankTime
                  << ",\"point_time\":" << p.pointTime << "}";
            }
            o << "]}" << std::endl;
        }
    }

    sum.ok = true;
    sum.numPatches = (int)patches.size();
    sum.flankRegions = flankRegions;
    sum.flankCut = flankCut;
    sum.flankOverhead = flankOverhead;
    sum.flankTotal = flankTotal;
    sum.pointCut = pointCut;
    sum.pointTotal = pointTotal;
    sum.speedup = speedup;
    sum.totalArea = totalArea;
    sum.originalArea = originalArea > 0.0 ? originalArea : 0.0;
    return sum;
}

} // namespace simple
