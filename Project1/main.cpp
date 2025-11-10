#include <glad/glad.h>
#include <glad/glad.c> // NOWE — musi byæ PRZED GLFW
#include <GLFW/glfw3.h>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <thread>           // std::thread, hardware_concurrency
#include <algorithm>        // std::sort, std::unique, std::min/std::max
#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
//#define DR_WAV_IMPLEMENTATION
//#include <dr_libs-master/dr_wav.h>
#include <iomanip>
#include <sstream>
#include <cctype>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_opengl3_loader.h"
#include <windows.h>
#define M_PI 3.1415

int test = 0;
bool simulation_running = false;
std::string file = "new_signal.csv";
std::string write_file = "new_signal_out.csv";
int windows_number = 50;

int ktore_odbicie = 0;
float prev_T2_time = 0;
float T2 = 0;

// gdzieœ przy sta³ych globalnych
static constexpr float R0_ATTEN = 0.05f; // [m] – near-field cap (np. 5 cm)

// — globals:
static uint32_t gFrameId = 0;
static float    gAvgEdgeLen = 0.0f;   // policzymy z siatki
static int      gBlindRings = 5;      // „3 trójk¹ty”
static float    gBlindRadius = 0.0f;  // = gBlindRings * gAvgEdgeLen

static constexpr float SOUND_V = 1440.0f / 1.0f; // C++
bool serio_first = true;
bool rewind_punkt = false;

float radius = 0.5f;
//float mic_radius = 0.2f;
const float H_ANGLE = M_PI / 180 * 72; // 72 stopni w radianach
const float V_ANGLE = atanf(1.0f / 2); // Kat wierzcholka

static size_t gWinIdx = 0;          // który 5 ms segment aktualnie nadajemy

//WAZNE DANE
float dt = 0.00001f;
float mic_radius = 0.5f;
float src_radius = 0.5f;
float time_passed = 0.0f;
float window_ms = 5.0f;
// --- NOWE: surowe dane z CSV ---
std::vector<double> tSec;       // czas w sekundach (po przeskalowaniu)
std::vector<float>  y;          // wartoœci próbek (2. kolumna)

// do podzialu sfery na klatki
// ==== PERSISTENT STATE dla rafinowania „po kawa³ku” ====
static size_t gRefineCursor = 0; // gdzie skoñczyliœmy ostatnio
static std::unordered_map<uint64_t, int> gEdgeMidCache; // edge -> midpoint

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

bool firstMouse = true;
bool mousePressed = false;
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
float yaw = -90.0f;
float pitch = 0.0f;
float fov = 45.0f;
float deltaTime = 0.0f;
float lastFrame = 0.0f;

static bool gMeshDirty = true;           // trzeba odbudowaæ bufory (zmiana triangles/nodes count)
struct MeshGL {
    GLuint  vbo = 0;
    GLuint  ibo = 0;
    GLsizei indexCount = 0;
    bool    dynamic = false;
};

struct Triangle {
    int indices[3];
};
std::vector<Triangle> triangles;
std::vector<Triangle> microphone;
std::vector<Triangle> source_tri;

struct WindowPacket {
    float tEmit = 0.0f;              // czas startu okna (sek)
    float amplitude = 0.0f;          // np. max |y| z okna (do progowania/prune)
    std::vector<float> times;        // czasy próbek RELATYWNE do tEmit [s], 0..5ms
    std::vector<float> values;       // wartoœci próbek (surowe)
};
static std::vector<WindowPacket> gWinPackets;

// trzy zestawy: fala, mikrofon, Ÿród³o
static MeshGL gWaveGL, gMicGL, gSrcGL;

// do FPS
int frameCount = 0;
double previousTime = 0.0;
double fps = 0.0;

bool first = true;
bool doKill = false;


struct MicSample { float t; float value; };   // czas (s), wartoœæ próbki
static std::vector<MicSample> gMicEvents;
std::vector<float> winMean;     // œrednie z okien X ms

struct node {
    glm::vec3 position = { 0,0,0 };
    glm::vec3 velocity = { 0,0,0 };    // kierunek propagacji * C_SOUND
    float     energy = 0.0f;
    uint8_t   bounces = 0; // w ostatecznej wersji mozna usunac chyba
    float     tEmit = 0.0f;       // czas emisji (sim-time)
    //int       winId = -1;   // <--- NOWE: identyfikator pakietu 5 ms
    // do Dopplera
    glm::vec3 srcVel = { 0,0,0 };    // prêdkoœæ Ÿród³a przy emisji
    glm::vec3 nEmit = { 0,0,0 };    // kierunek promienia przy emisji (unit)

    int       seedId = -1;            // indeks wierzcho³ka siatki (opcjonalnie)
    //uint32_t  suppressStamp = 0;      // znacznik „nie zapisuj do CSV” w tym kroku
    float     suppressUntilT = -1e30f; // NOWE: nie logowaæ, jeœli T_hit <= suppressUntilT

    float path = 0.0f;
    float doppler = 1.0f;
};

std::vector<node> nodes;
std::vector<node> mic_nodes;
std::vector<node> src_nodes;

struct source {
    float src_x = 0.0f;
    float src_y = 0.0f;
    float src_z = 0.0f;
    //glm::vec3 starting_point = glm::vec3(src_x, src_y, src_z);
    glm::vec3 velocity = glm::vec3(-00.0f, 0, 0);
    glm::vec3 rewind_point = glm::vec3(src_x, src_y, src_z);
    glm::vec3 rewind_vel = velocity;
};
source Source;

struct Cuboid_dimensions {
    float width = 0.0f;
    float height = 0.0f;
    float depth = 0.0f;
    float x_offset = 0.0f;
    float y_offset = 0.0f;
    float z_offset = 0.0f;
};
Cuboid_dimensions Cube{
    100.0f, 100.0f, 100.0f,  // width, height, depth
    0.0f,   0.0f,    0.0f        // x/y/z offset
};

Cuboid_dimensions Obstacle{
    0.3f, 0.7f, 0.7f,       //width, height, depth
    50.0f, 1.0f, 50.0f      // x/y/z offset
};

//MIKROFON
struct Micophone {
    float mic_x = 1.5f;
    float mic_y = 0.0f;
    float mic_z = 0.0f;
    //glm::vec3 starting_point = glm::vec3(mic_x, mic_y, mic_z);
    glm::vec3 mic_velocity = glm::vec3(100.0f, 000.0f, 0.0f);
    glm::vec3 rewind_point = glm::vec3(mic_x, mic_y, mic_z);
    glm::vec3 rewind_vel = mic_velocity;
};
Micophone Mic;




//#include <cstdlib> // wymagane dla exit()
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void renderScene();
//void drawCuboid(float width, float height, float depth);
void drawCuboidTransparentSorted(struct Cuboid_dimensions temp_Cube);
int printOversizedTriangles(float maxArea);
// Ustawia wêz³y w pozycji Ÿród³a i nadaje im prêdkoœci/kierunki startowe.
void killAllNodes();
void buildBuffersFor(const std::vector<node>& verts, const std::vector<Triangle>& tris, MeshGL& m, bool dynamic = true);
void updatePositionsFor(const std::vector<node>& verts, MeshGL& m);
void drawMesh(const MeshGL& m, const glm::vec3& offset = glm::vec3(0), const glm::vec4& fill = glm::vec4(1.0f, 0.5f, 0.0f, 1.0f), const glm::vec4& wire = glm::vec4(0.0f, 0.0f, 0.0f, 0.6f), float wireWidth = 0.2f);
static void buildIcosphereNodesTris(float radius, int subdiv, std::vector<node>& out_nodes, std::vector<Triangle>& out_tris);
//Zapis do pliku po symulacji
//int removeSlowNodes(float minSpeed);

static inline float q5(float t) {
    // zaokr¹glenie do najbli¿szego kroku 1e-5 s (5 miejsc po przecinku)
    // double w œrodku zmniejsza b³êdy kumulacji
    return (float)(std::llround((double)t * 100000.0) / 100000.0);
}

static float computeAvgEdgeLen(const std::vector<glm::vec3>& V,
    const std::vector<glm::ivec3>& F) {
    double sum = 0; size_t cnt = 0;
    for (auto& t : F) {
        const auto& a = V[t.x], & b = V[t.y], & c = V[t.z];
        sum += glm::length(a - b); ++cnt;
        sum += glm::length(b - c); ++cnt;
        sum += glm::length(c - a); ++cnt;
    }
    return cnt ? float(sum / cnt) : 0.0f;
}

// Autodetekcja separatora: wybierz ';' jeœli jest, inaczej ','
static inline char detect_sep(const std::string& line) {
    return (line.find(';') != std::string::npos) ? ';' : ',';
}

inline void logMicHit(float t_dop, float value) {
    if (t_dop < 0.0f) t_dop = 0.0f;          // na wszelki wypadek
    // std::lock_guard<std::mutex> lk(gMicMtx);
    gMicEvents.push_back({ t_dop, value });
}

inline void resetMicEvents() { gMicEvents.clear(); }


bool writeMicCsv(const std::string& path = write_file) {
    if (gMicEvents.empty()) return false;

    // 1) posortuj po czasie
    std::sort(gMicEvents.begin(), gMicEvents.end(),
        [](const MicSample& a, const MicSample& b) { return a.t < b.t; });

    // 2) scal duplikaty czasów (dok³adnie równe) przez sumowanie wartoœci
    std::vector<MicSample> merged;
    merged.reserve(gMicEvents.size());
    for (const auto& s : gMicEvents) {
        if (!merged.empty() && merged.back().t == s.t) {
            merged.back().value += s.value;
        }
        else {
            merged.push_back(s);
        }
    }

    // 3) zapisz do CSV
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    if (!f) return false;

    f.setf(std::ios::fixed);
    f << "time,value\n";
    f << std::setprecision(9);
    for (const auto& s : merged) {
        f << s.t << ",";
        f << std::setprecision(6) << s.value << "\n";
        f << std::setprecision(9); // wróæ do prec dla czasu
    }
    return true;
}

// Parsuj CSV: kol.1 = czas, kol.2 = wartoœæ.
// prosta wersja: bez sortowania, autodetekcji separatora itp.
bool loadCsvSimple(const std::string& path,
    char sep = ',',          // ustaw ';' jeœli tak masz
    bool hasHeader = true,
    double timeScale = 1.0)  // 1.0 gdy time_s w sekundach; 0.001 gdy w ms
{
    tSec.clear(); y.clear(); winMean.clear();

    std::ifstream f(path);
    if (!f) return false;

    std::string line;
    if (hasHeader) std::getline(f, line); // pomiñ nag³ówek
    std::cout << line << std::endl;

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        size_t pos = line.find(sep);
        if (pos == std::string::npos) continue;

        double t = std::stod(line.substr(0, pos)) * timeScale; // -> sekundy
        double v = std::stod(line.substr(pos + 1));
        tSec.push_back(t);
        y.push_back((float)v);
    }
    if (tSec.empty()) return false;

    // Sta³e 5 ms
    const double win_s = 0.005;          // 5 ms w sekundach
    const double t0 = tSec.front();
    const double tEnd = tSec.back();
    const size_t nWin = (size_t)std::floor((tEnd - t0) / win_s) + 1;

    // Przygotuj pakiety i (dla zgodnoœci) winMean jako "amplitudê"
    gWinPackets.clear();
    gWinPackets.resize(nWin);
    winMean.assign(nWin, 0.0f);
    std::vector<float> maxAbs(nWin, 0.0f);

    // Zainicjuj tEmit dla ka¿dego okna
    for (size_t k = 0; k < nWin; ++k) {
        gWinPackets[k].tEmit = float(t0 + k * win_s);
    }

    // Rozdziel próbki do okien
    for (size_t i = 0; i < tSec.size(); ++i) {
        size_t k = (size_t)((tSec[i] - t0) / win_s);
        if (k >= nWin) k = nWin - 1;

        float tRel = float(tSec[i] - (t0 + k * win_s));  // czas wzgl. startu okna
        gWinPackets[k].times.push_back(tRel);
        gWinPackets[k].values.push_back(y[i]);

        float a = std::fabs(y[i]);
        if (a > maxAbs[k]) maxAbs[k] = a;
    }

    // Ustal amplitude i dla zgodnoœci wpisz j¹ te¿ do winMean
    for (size_t k = 0; k < nWin; ++k) {
        gWinPackets[k].amplitude = maxAbs[k]; // definicja amplitudy: max|y|
        winMean[k] = maxAbs[k];               // ¿eby getAtTime nadal dawa³ sensowny próg
    }
    return true;
}

// Œrednia okna odpowiadaj¹cego czasowi t [s]
float getAtTime(double t_sec) {
    if (!gWinPackets.empty()) {
        const double win_s = 0.005; // 5 ms
        size_t k = (size_t)std::floor((t_sec + 1e-9) / win_s);
        if (k >= gWinPackets.size()) return 0.0f;
        return gWinPackets[k].amplitude; // amplituda okna
    }
    // Fallback do starego
    if (winMean.empty()) return 0.0f;
    const double idx = (t_sec * 1000.0 + 1e-3) / double(window_ms);
    size_t k = (size_t)std::floor(idx);
    if (k >= winMean.size()) return 0.0f;
    return winMean[k];
}


static bool beginNextWindow() {
    if (gWinIdx >= gWinPackets.size()) return false;
    gWinIdx++;
    first = true;   // wypuszczamy now¹ falê
    return true;
}

int pruneSlowNodes(float minSpeed);

static inline bool touchesMicrophone(const glm::vec3& p);
// unikalny klucz krawêdzi (mniejszy indeks najpierw)
static inline uint64_t edge_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (uint64_t(a) << 32) | uint32_t(b);
}

// rzut na sferê (promieñ r)
static inline glm::vec3 normalize_to_radius(const glm::vec3& p, float r) {
    float len = glm::length(p);
    if (len == 0.0f) return glm::vec3(r, 0, 0);
    return (p / len) * r;
}

// zwraca indeks œrodkowego wierzcho³ka (z cache'owaniem, by unikn¹æ duplikatów)
static int midpoint_index(int i0, int i1,
    std::vector<glm::vec3>& verts,
    std::unordered_map<uint64_t, int>& cache,
    float radius)
{
    uint64_t key = edge_key(i0, i1);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    glm::vec3 m = 0.5f * (verts[i0] + verts[i1]);
    m = normalize_to_radius(m, radius);       // PROJEKCJA NA SFERÊ
    int idx = (int)verts.size();
    verts.push_back(m);
    cache.emplace(key, idx);
    return idx;
}

void updatePhysics(float dt, struct Cuboid_dimensions Pool, struct Cuboid_dimensions temp_Obstacle)
{
    const float Pool_halfW = 0.5f * Pool.width;
    const float Pool_halfH = 0.5f * Pool.height;
    const float Pool_halfD = 0.5f * Pool.depth;

    const float e = 0.8;   // wsp. sprê¿ystoœci
    const float eps = 0.001f; // minimalne odsuniêcie od œciany

    const float micR = mic_radius;
    ++gFrameId;
    // Pomocnik do odbicia w 1D (szybki, inline'owalny)
    auto bounce1D = [&](float& pos, float& vel, float minb, float maxb) -> bool
        {
            if (pos < minb) {
                pos = minb + eps;    // wyci¹gnij z penetracji
                vel = -vel * e;      // odbij (uwzglêdnij stratê)
                return true;
            }
            else if (pos > maxb) {
                pos = maxb - eps;
                vel = -vel * e;
                return true;
            }
            return false;
        };

    auto bounceObstacleMic = [&](struct Micophone& temp_Mic, struct Cuboid_dimensions temp_Obstacle)
        {
            const float temp_Obstacle_halfW = 0.5f * temp_Obstacle.width;
            const float temp_Obstacle_halfH = 0.5f * temp_Obstacle.height;
            const float temp_Obstacle_halfD = 0.5f * temp_Obstacle.depth;

            if (temp_Obstacle_halfW + temp_Obstacle.x_offset >= 0)
            {
                if (temp_Mic.mic_x > -temp_Obstacle_halfW + temp_Obstacle.x_offset && temp_Mic.mic_x < temp_Obstacle_halfW + temp_Obstacle.x_offset)
                {
                    temp_Mic.mic_velocity.x *= -1;
                }
            }
            else
            {
                if (temp_Mic.mic_x < -temp_Obstacle_halfW + temp_Obstacle.x_offset && temp_Mic.mic_x > temp_Obstacle_halfW + temp_Obstacle.x_offset)
                {
                    temp_Mic.mic_velocity.x *= -1;
                }
            }

            if (temp_Obstacle_halfH + temp_Obstacle.y_offset >= 0)
            {
                if (temp_Mic.mic_y > -temp_Obstacle_halfH + temp_Obstacle.y_offset && temp_Mic.mic_y < temp_Obstacle_halfH + temp_Obstacle.y_offset)
                {
                    temp_Mic.mic_velocity.y *= -1;
                }
            }
            else
            {
                if (temp_Mic.mic_y < -temp_Obstacle_halfH + temp_Obstacle.y_offset && temp_Mic.mic_y > temp_Obstacle_halfH + temp_Obstacle.y_offset)
                {
                    temp_Mic.mic_velocity.y *= -1;
                }
            }

            if (temp_Obstacle_halfD + temp_Obstacle.z_offset >= 0)
            {
                if (temp_Mic.mic_z > -temp_Obstacle_halfD + temp_Obstacle.z_offset && temp_Mic.mic_z < temp_Obstacle_halfD + temp_Obstacle.z_offset)
                {
                    temp_Mic.mic_velocity.z *= -1;
                }
            }
            else
            {
                if (temp_Mic.mic_z < -temp_Obstacle_halfD + temp_Obstacle.z_offset && temp_Mic.mic_z > temp_Obstacle_halfD + temp_Obstacle.z_offset)
                {
                    temp_Mic.mic_velocity.z *= -1;
                }
            }

        };

    auto bounceObstacleWave = [&](node& temp_node, struct Cuboid_dimensions temp_Obstacle)
        {
            const float temp_Obstacle_halfW = 0.5f * temp_Obstacle.width;
            const float temp_Obstacle_halfH = 0.5f * temp_Obstacle.height;
            const float temp_Obstacle_halfD = 0.5f * temp_Obstacle.depth;

            if ((temp_node.position.x > -temp_Obstacle_halfW + temp_Obstacle.x_offset && temp_node.position.y > -temp_Obstacle_halfH + temp_Obstacle.y_offset
                && temp_node.position.z > -temp_Obstacle_halfD + temp_Obstacle.z_offset)
                &&
                (temp_node.position.x < temp_Obstacle_halfW + temp_Obstacle.x_offset && temp_node.position.y < temp_Obstacle_halfH + temp_Obstacle.y_offset
                    && temp_node.position.z < temp_Obstacle_halfD + temp_Obstacle.z_offset))
            {
                //std::cout << "KOLIZJA Z FALA" << std::endl;
                if (temp_node.position.x - 5 * eps < -temp_Obstacle_halfW + temp_Obstacle.x_offset || temp_node.position.x + 5 * eps > temp_Obstacle_halfW + temp_Obstacle.x_offset)
                {
                    temp_node.velocity.x *= -1;
                }
                else if (temp_node.position.y - 5 * eps < -temp_Obstacle_halfH + temp_Obstacle.y_offset || temp_node.position.y + 5 * eps > temp_Obstacle_halfH + temp_Obstacle.y_offset)
                {
                    temp_node.velocity.y *= -1;
                }
                else if (temp_node.position.z - 5 * eps < -temp_Obstacle_halfD + temp_Obstacle.z_offset || temp_node.position.z + 5 * eps > temp_Obstacle_halfD + temp_Obstacle.z_offset)
                {
                    temp_node.velocity.z *= -1;
                }
            }


        };


    // --- 1) Integracja i odbicie mikrofonu (poza pêtl¹ równoleg³¹) ---
    Mic.mic_x += Mic.mic_velocity.x * dt;
    Mic.mic_y += Mic.mic_velocity.y * dt;
    Mic.mic_z += Mic.mic_velocity.z * dt;

    Source.src_x += Source.velocity.x * dt;
    Source.src_y += Source.velocity.y * dt;
    Source.src_z += Source.velocity.z * dt;

    // Odbicia mikrofonu od œcian basenu, z uwzglêdnieniem promienia
    bounce1D(Mic.mic_x, Mic.mic_velocity.x, -Pool_halfW + Pool.x_offset + micR, Pool_halfW + Pool.x_offset - micR);
    bounce1D(Mic.mic_y, Mic.mic_velocity.y, -Pool_halfH + Pool.y_offset + micR, Pool_halfH + Pool.y_offset - micR);
    bounce1D(Mic.mic_z, Mic.mic_velocity.z, -Pool_halfD + Pool.z_offset + micR, Pool_halfD + Pool.z_offset - micR);

    //odbicia od przeszkody (MIKROFON)
    bounceObstacleMic(Mic, temp_Obstacle);

    //doKill = (glfwGetTime() >= 8.0);
    // --- 2) Integracja i odbicia punktów siatki ---
//#pragma omp parallel for schedule(static)
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        auto& p = nodes[i].position;
        auto& v = nodes[i].velocity;
        auto& energy = nodes[i].energy;

        //nodes[i].energy *= e;
        //energia_test = nodes[i].energy;

        // nowe pozycje
        p += v * dt;
        // d³ugoœæ przebytego odcinka w tej klatce:
        float stepLen = glm::length(v) * dt;           // (poprawka GLM!)
        float r_prev = std::max(nodes[i].path, R0_ATTEN);
        nodes[i].path += stepLen * 0.01;
        float r_now = std::max(nodes[i].path, R0_ATTEN);

        // energia ~ 1/r  (utrzymuje proporcjê niezale¿nie od kroku dt)
        energy *= (r_prev / r_now);
        //energy = gWinPackets[gWinIdx].amplitude / nodes[i].path;
        //energy = energy * 0.999999;


        // Odbicia w XYZ
        // --- w pêtli po node'ach ---
        bool bouncedAny = false;

        if (bounce1D(p.x, v.x, -Pool_halfW + Pool.x_offset, Pool_halfW + Pool.x_offset))
        {
            bouncedAny |= true;
            energy = energy * 0.5;
        }
        if (bounce1D(p.y, v.y, -Pool_halfH + Pool.y_offset, Pool_halfH + Pool.y_offset))
        {
            bouncedAny |= true;
            energy = energy * 0.5;
        }
        if (bounce1D(p.z, v.z, -Pool_halfD + Pool.z_offset, Pool_halfD + Pool.z_offset))
        {
            bouncedAny |= true;
            energy = energy * 0.5;
        }

        //bouncedAny |= bounce1D(p.x, v.x, -Pool_halfW + Pool.x_offset, Pool_halfW + Pool.x_offset);
        //bouncedAny |= bounce1D(p.y, v.y, -Pool_halfH + Pool.y_offset, Pool_halfH + Pool.y_offset);
        //bouncedAny |= bounce1D(p.z, v.z, -Pool_halfD + Pool.z_offset, Pool_halfD + Pool.z_offset);

        if (bouncedAny) {
            nodes[i].bounces++;      // +1 za "jakiekolwiek" odbicie w tym kroku
        }

        //------ODBICIA ZRODLA OD SCIAN------
        bounce1D(Source.src_x, Source.velocity.x, -Pool_halfW + Pool.x_offset, Pool_halfW + Pool.x_offset);
        bounce1D(Source.src_y, Source.velocity.y, -Pool_halfH + Pool.y_offset, Pool_halfH + Pool.y_offset);
        bounce1D(Source.src_z, Source.velocity.z, -Pool_halfD + Pool.z_offset, Pool_halfD + Pool.z_offset);


        //odbicia od przeszkody (FALA)
        //--------MOZNA ZAKOMENTOWAC--------
        bounceObstacleWave(nodes[i], temp_Obstacle);

    }
    //if (doKill) doKill = false;
    //dodaj czas
    time_passed += dt;
    //if (time_passed > 0.0008)
    //{
    //    int zero = 0;
    //}

    if (time_passed / dt >= (window_ms / 1000.0f) / dt && rewind_punkt)
    {
        Mic.rewind_point = glm::vec3(Mic.mic_x, Mic.mic_y, Mic.mic_z);
        Mic.rewind_vel = Mic.mic_velocity;
        Source.rewind_point = glm::vec3(Source.src_x, Source.src_y, Source.src_z);
        Source.rewind_vel = Source.velocity;
        rewind_punkt = false;
    }

    //rewind_punkt = false;

}

float calculateTriangleArea(int a, int b, int c) {
    glm::vec3 ab = nodes[b].position - nodes[a].position;
    glm::vec3 ac = nodes[c].position - nodes[a].position;
    glm::vec3 cross = glm::cross(ab, ac);
    return 0.5f * glm::length(cross); //funkcja prawdziwego pola
    //return glm::dot(cross, cross); // pole bez pierwiastka czyli (2*area)^2
}

int addMidpoint(int a, int b) {
    node midpoint;
    midpoint.position = (nodes[a].position + nodes[b].position) * 0.5f;
    midpoint.velocity = (nodes[a].velocity + nodes[b].velocity) * 0.5f;
    midpoint.energy = (nodes[a].energy + nodes[b].energy) * 0.5f;
    midpoint.tEmit = (nodes[a].tEmit + nodes[b].tEmit) * 0.5f;
    midpoint.srcVel = (nodes[a].srcVel + nodes[b].srcVel) * 0.5f;
    midpoint.bounces = nodes[a].bounces;
    midpoint.suppressUntilT = (nodes[a].suppressUntilT + nodes[b].suppressUntilT) * 0.5f;
    midpoint.path = (nodes[a].path + nodes[b].path) * 0.5f;
    //midpoint.velocity = (nodes[a].velocity + nodes[b].velocity) * 0.5f;

    nodes.push_back(midpoint);
    return (int)(nodes.size() - 1);
}

// ====== POMOCNICZE ======
struct Edge { int a, b; }; // zawsze a<b
static inline Edge make_edge(int u, int v) { if (u > v) std::swap(u, v); return { u,v }; }
static inline bool edge_less(const Edge& x, const Edge& y) {
    return (x.a < y.a) || (x.a == y.a && x.b < y.b);
}
static inline bool edge_eq(const Edge& x, const Edge& y) {
    return x.a == y.a && x.b == y.b;
}

static int midpoint_index_nodes_cached(int i0, int i1) {
    uint64_t key = edge_key(i0, i1);           // masz ju¿ edge_key(i,j)
    auto it = gEdgeMidCache.find(key);
    if (it != gEdgeMidCache.end()) return it->second;
    int idx = addMidpoint(i0, i1);      // Twój kod midpointu (pozycja/prêdkoœæ)
    gEdgeMidCache.emplace(key, idx);
    return idx;
}

// Zwraca indeks midpointu z cache; jeœli nie istnieje — tworzy przez addMidpoint(...)
static int midpoint_index_nodes(int i0, int i1,
    std::unordered_map<uint64_t, int>& cache)
{
    uint64_t key = edge_key(i0, i1);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    int idx = addMidpoint(i0, i1); // œrednia pozycji i prêdkoœci
    cache.emplace(key, idx);
    return idx;
}
// ====== MULTI-THREAD CHUNKED REFINEMENT ======
bool refineIcosahedron_chunked_mt(float maxArea,
    size_t triBudget,          // ile trójk¹tów obrabiamy w tej klatce
    int   threadCount = std::thread::hardware_concurrency())
{
    if (triangles.empty() || triBudget == 0) return false;
    if (threadCount <= 0) threadCount = 1;

    // jeœli siatka by³a resetowana
    if (gRefineCursor > triangles.size()) {
        gRefineCursor = 0;
        gEdgeMidCache.clear();
    }

    // próg bez sqrt: |cross|^2 > (2*maxArea)^2 zamiast porownywac pola porownywamy kwadraty 
    const float area2_threshold = 4.0f * maxArea * maxArea;

    // sta³e kolizyjne
    const float cuboidHalfWidth = Cube.width * 0.5f;
    const float cuboidHalfHeight = Cube.height * 0.5f;
    const float cuboidHalfDepth = Cube.depth * 0.5f;
    const float safetyMargin = 0.2f;

    // pracujemy wy³¹cznie na „starym” prefiksie tablicy trójk¹tów
    const size_t size0 = triangles.size();
    size_t start = (gRefineCursor < size0 ? gRefineCursor : 0);
    size_t end = std::min(start + triBudget, size0);
    const size_t N = (end > start ? end - start : 0);
    if (N == 0) return false;

    threadCount = std::min<int>(threadCount, (int)N);

    // --- FAZA 1: RÓWNOLEGLE zbierz „do podzia³u” + krawêdzie ---
    std::vector<std::vector<int>>   tl_split(threadCount);
    std::vector<std::vector<Edge>>  tl_edges(threadCount);

    auto worker_scan = [&](int tid) {
        size_t chunk = (N + threadCount - 1) / threadCount;
        size_t i0 = start + tid * chunk;
        size_t i1 = std::min(end, i0 + chunk);
        auto& split = tl_split[tid];
        auto& edges = tl_edges[tid];
        split.reserve((i1 > i0) ? (i1 - i0) / 2 : 0);
        edges.reserve((i1 > i0) ? 3 * (i1 - i0) / 2 : 0);

        for (size_t i = i0; i < i1; ++i) {
            Triangle t = triangles[i]; // kopiuj (nie referencja; póŸniej bêdziemy pushowaæ do triangles)
            int a = t.indices[0], b = t.indices[1], c = t.indices[2];

            const glm::vec3& pa = nodes[a].position;
            const glm::vec3& pb = nodes[b].position;
            const glm::vec3& pc = nodes[c].position;

            const bool nearWall =
                (std::fabs(pa.x) > cuboidHalfWidth - safetyMargin) ||
                (std::fabs(pa.y) > cuboidHalfHeight - safetyMargin) ||
                (std::fabs(pa.z) > cuboidHalfDepth - safetyMargin) ||
                (std::fabs(pb.x) > cuboidHalfWidth - safetyMargin) ||
                (std::fabs(pb.y) > cuboidHalfHeight - safetyMargin) ||
                (std::fabs(pb.z) > cuboidHalfDepth - safetyMargin) ||
                (std::fabs(pc.x) > cuboidHalfWidth - safetyMargin) ||
                (std::fabs(pc.y) > cuboidHalfHeight - safetyMargin) ||
                (std::fabs(pc.z) > cuboidHalfDepth - safetyMargin);

            if (nearWall) continue;

            const glm::vec3 ab = pb - pa;
            const glm::vec3 ac = pc - pa;
            const glm::vec3 cr = glm::cross(ab, ac);
            const float     cr2 = glm::dot(cr, cr);

            if (cr2 > area2_threshold) {
                split.push_back((int)i);
                edges.push_back(make_edge(a, b));
                edges.push_back(make_edge(b, c));
                edges.push_back(make_edge(c, a));
            }
        }
        };

    std::vector<std::thread> threads;
    threads.reserve(threadCount);
    for (int t = 0; t < threadCount; ++t) threads.emplace_back(worker_scan, t);
    for (auto& th : threads) th.join();
    threads.clear();

    // scalenie wyników fazy 1
    std::vector<int> splitIdx;
    std::vector<Edge> edges;
    {
        size_t totalSplit = 0, totalEdges = 0;
        for (int t = 0; t < threadCount; ++t) { totalSplit += tl_split[t].size(); totalEdges += tl_edges[t].size(); }
        if (totalSplit == 0) { // nic do roboty
            gRefineCursor = end < size0 ? end : 0;
            return false;
        }
        splitIdx.reserve(totalSplit);
        edges.reserve(totalEdges);
        for (int t = 0; t < threadCount; ++t) {
            splitIdx.insert(splitIdx.end(), tl_split[t].begin(), tl_split[t].end());
            edges.insert(edges.end(), tl_edges[t].begin(), tl_edges[t].end());
            tl_split[t].clear(); tl_edges[t].clear();
        }
    }

    // unikalizacja krawêdzi
    std::sort(edges.begin(), edges.end(), edge_less);
    edges.erase(std::unique(edges.begin(), edges.end(), edge_eq), edges.end());

    // --- FAZA 2: SEKWENCYJNIE utwórz midpointy (u¿yj globalnego cache) ---
    std::vector<int> edgeMidIdx(edges.size());
    nodes.reserve(nodes.size() + edges.size());          // minimalizacja realokacji
    for (size_t i = 0; i < edges.size(); ++i) {
        edgeMidIdx[i] = midpoint_index_nodes_cached(edges[i].a, edges[i].b);
    }
    //auto edge_lookup = [&](int u, int v)->int {
    //    if (u > v) std::swap(u, v);
    //    Edge key{ u,v };
    //    auto it = std::lower_bound(edges.begin(), edges.end(), key, edge_less);
    //    // przy poprawnym zbiorze powinien istnieæ
    //    return edgeMidIdx[size_t(it - edges.begin())];
    //    };

    auto edge_lookup = [&](int u, int v) -> int {
        if (u > v) std::swap(u, v);
        const uint64_t key = edge_key(u, v);
        auto it = gEdgeMidCache.find(key);
        return it->second;
        };

    // --- FAZA 3: RÓWNOLEGLE zbuduj nowe trójk¹ty do buforów per-w¹tek ---
    std::vector<std::vector<std::pair<int, Triangle>>> tl_replace(threadCount);
    std::vector<std::vector<Triangle>>                tl_append(threadCount);

    auto worker_build = [&](int tid) {
        size_t M = splitIdx.size();
        size_t chunk = (M + threadCount - 1) / threadCount;
        size_t i0 = tid * chunk;
        size_t i1 = std::min(M, i0 + chunk);
        auto& rep = tl_replace[tid];
        auto& app = tl_append[tid];
        rep.reserve((i1 > i0) ? (i1 - i0) : 0);
        app.reserve((i1 > i0) ? 3 * (i1 - i0) : 0);

        for (size_t k = i0; k < i1; ++k) {
            int ti = splitIdx[k];
            Triangle t = triangles[ti];
            int a = t.indices[0], b = t.indices[1], c = t.indices[2];

            int ab = edge_lookup(a, b);
            int bc = edge_lookup(b, c);
            int ca = edge_lookup(c, a);

            // zamieniamy bie¿¹cy
            rep.emplace_back(ti, Triangle{ {a,  ab, ca} });
            // 3 nowe
            app.push_back({ {b,  bc, ab} });
            app.push_back({ {c,  ca, bc} });
            app.push_back({ {ab, bc, ca} });
        }
        };

    for (int t = 0; t < threadCount; ++t) threads.emplace_back(worker_build, t);
    for (auto& th : threads) th.join();

    // --- FAZA 4: JEDNYM KROKIEM zapisz wynik do 'triangles' ---
    // (po fazie 3 nie czytamy ju¿ 'triangles', wiêc mo¿emy rezerwowaæ i pushowaæ)
    size_t toAppend = 0;
    for (int t = 0; t < threadCount; ++t) toAppend += tl_append[t].size();
    triangles.reserve(triangles.size() + toAppend);

    // podmiany
    for (int t = 0; t < threadCount; ++t) {
        for (auto& p : tl_replace[t]) {
            triangles[p.first] = p.second;
        }
    }
    // dopisywanie 3 nowych
    for (int t = 0; t < threadCount; ++t) {
        triangles.insert(triangles.end(), tl_append[t].begin(), tl_append[t].end());
    }

    // przesuwamy kursor (zawsze wzglêdem dawnego size0)
    gRefineCursor = end;
    if (gRefineCursor >= size0) gRefineCursor = 0;

    return true;
}



void calculateFPS() {
    double currentTime = glfwGetTime();
    frameCount++;

    if (currentTime - previousTime >= 1.0) {
        fps = frameCount / (currentTime - previousTime);
        std::cout << "FPS: " << fps << std::endl;

        frameCount = 0;
        previousTime = currentTime;
        printOversizedTriangles(0.05f);
    }
}


void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if (action == GLFW_PRESS) {
            mousePressed = true;
        }
        else if (action == GLFW_RELEASE) {
            mousePressed = false;
            firstMouse = true;
        }
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    // Zmiana odleglosci kamery za pomoca scrolla
    cameraPos += cameraFront * static_cast<float>(yoffset) * 0.5f;
}

static void glfwErrorCallback(int code, const char* desc) { fprintf(stderr, "GLFW[%d]: %s\n", code, desc); }

//---------------------
//--------MAIN---------
//---------------------
int main()
{


    glfwSetErrorCallback(glfwErrorCallback);
    if (!glfwInit())
        return -1;

    bool pokazano_dane = false; //testowe, do usuniecia

    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Camera Scene", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);



    // Wczytaj plik WAV (podmieñ œcie¿kê na swoj¹)
    if (!loadCsvSimple(file, ',', /*hasHeader=*/true)) {
        std::cerr << "CSV: nie wczytano ¿adnych danych\n";
    }
    else {
        std::cout << "CSV: wczytano " << tSec.size()
            << " wierszy" << std::endl;

    }
    //beginNextWindow(); // uruchamiamy pierwsze okno 5 ms
    //gAudio.window_ms = 5.0f;   // trzymamy 5 ms

    // to pod tym dodane do VBO 
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "GLAD fail\n"; return -1;
    }
    std::cout << "GL " << GLVersion.major << "." << GLVersion.minor << "\n";
    //glEnable(GL_DEPTH_TEST);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);


    //Rzeczy od ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext(); //inicjalizacja ImGui     
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark(); //ciemny motyw imgui
    ImGui_ImplGlfw_InitForOpenGL(window, true); //setup
    ImGui_ImplOpenGL3_Init("#version 330");



    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);

        glClearColor(0.5f, 0.5f, 0.5f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 projection = glm::perspective(glm::radians(fov), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);

        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(&projection[0][0]);

        glMatrixMode(GL_MODELVIEW);
        glLoadMatrixf(&view[0][0]);


        // Oblicz FPS
        //calculateFPS();
        //if (simulation_running)
        //{
        renderScene();
        //}

        glfwPollEvents();

        //---IMGUI---
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        ImGui::SetWindowSize(ImVec2(200, 200));

        ImGui::Begin("Start/Stop");
        if (ImGui::Button("Start"))
        {
            simulation_running = true;
        }
        if (ImGui::Button("Stop"))
        {
            simulation_running = false;
        }
        ImGui::End();

        ImGui::Begin("Parametry");
        ImGui::Text("Mikrofon");
        ImGui::InputFloat("Mic x", &Mic.mic_x, 0.1f, Mic.mic_x);
        ImGui::InputFloat("Mic y", &Mic.mic_y, 0.1f, Mic.mic_y);
        ImGui::InputFloat("Mic z", &Mic.mic_z, 0.1f, Mic.mic_z);
        ImGui::InputFloat("Mic_vel x", &Mic.mic_velocity.x, 5.0f, Mic.mic_velocity.x);
        ImGui::InputFloat("Mic_vel y", &Mic.mic_velocity.y, 5.0f, Mic.mic_velocity.y);
        ImGui::InputFloat("Mic_vel z", &Mic.mic_velocity.z, 5.0f, Mic.mic_velocity.z);
        ImGui::Text("Zrodlo");
        ImGui::InputFloat("Src x", &Source.src_x, 0.1f, Source.src_x);
        ImGui::InputFloat("Src y", &Source.src_y, 0.1f, Source.src_y);
        ImGui::InputFloat("Src z", &Source.src_z, 0.1f, Source.src_z);
        ImGui::InputFloat("Src_vel x", &Source.velocity.x, 5.0f, Source.velocity.x);
        ImGui::InputFloat("Src_vel y", &Source.velocity.y, 5.0f, Source.velocity.y);
        ImGui::InputFloat("Src_vel z", &Source.velocity.z, 5.0f, Source.velocity.z);
        ImGui::Text("Basen");
        ImGui::InputFloat("Cube width", &Cube.width, 2.0f, Cube.width);
        ImGui::InputFloat("Cube height", &Cube.height, 2.0f, Cube.height);
        ImGui::InputFloat("Cube depth", &Cube.depth, 2.0f, Cube.depth);
        ImGui::InputFloat("Cube x_offset", &Cube.x_offset, 2.0f, Cube.x_offset);
        ImGui::InputFloat("Cube y_offset", &Cube.y_offset, 2.0f, Cube.y_offset);
        ImGui::InputFloat("Cube z_offset", &Cube.z_offset, 2.0f, Cube.z_offset);
        ImGui::End();



        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        //---Koniec IMGUI----

        glfwSwapBuffers(window);


    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();


    glfwTerminate();
    return 0;
}

void processInput(GLFWwindow* window)
{
    float cameraSpeed = 2.5f * deltaTime;

    // Zmiana wysokosci kamery za pomoca klawiszy W i S
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos.y += cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos.y -= cameraSpeed;

    // Poruszanie kamera w plaszczyznie XZ za pomoca klawiszy A i D
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (!mousePressed)
        return;

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;

    if (pitch > 89.0f)
        pitch = 89.0f;
    if (pitch < -89.0f)
        pitch = -89.0f;

    glm::vec3 front;
    front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front.y = sin(glm::radians(pitch));
    front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraFront = glm::normalize(front);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}


void renderScene()
{
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //static std::vector<Triangle> triangles;
    if (first && !rewind_punkt) {
        doKill = false;
        //glfwSetTime(0.0);  // wyzeruj stoper
        time_passed = 0.0f;
        //microphone
        Mic.mic_x = Mic.rewind_point.x;
        Mic.mic_y = Mic.rewind_point.y;
        Mic.mic_z = Mic.rewind_point.z;
        Mic.mic_velocity = Mic.rewind_vel;
        //source
        Source.src_x = Source.rewind_point.x;
        Source.src_y = Source.rewind_point.y;
        Source.src_z = Source.rewind_point.z;
        Source.velocity = Source.rewind_vel;

        rewind_punkt = true;

        nodes.clear();
        triangles.clear();

        if (serio_first)
        {
            serio_first = false;
            // Zbuduj geometriê dla mikrofonu i Ÿród³a — z tym samym SUBDIV co fala,
           // lub jeœli chcesz, z osobnymi (np. SUBDIV_RIGID).
            const int SUBDIV_RIGID = 2; // mo¿esz daæ 0..2 niezale¿nie od fali

            mic_nodes.clear();
            microphone.clear();
            buildIcosphereNodesTris(mic_radius, SUBDIV_RIGID, mic_nodes, microphone);

            src_nodes.clear();
            source_tri.clear();
            buildIcosphereNodesTris(src_radius, SUBDIV_RIGID, src_nodes, source_tri);

            // Bufory GPU dla statycznych meshów (model-space):
            buildBuffersFor(mic_nodes, microphone, gMicGL, /*dynamic=*/false);
            buildBuffersFor(src_nodes, source_tri, gSrcGL, /*dynamic=*/false);
        }
        // --- Parametr gêstoœci: ka¿dy poziom ×4 liczba trójk¹tów ---
        constexpr int SUBDIV = 2; // 2 optymalnie, wiecej laguje
        // --- 12 wierzcho³ków  na sferze o promieniu 'radius' ---
        const float t = (1.0f + std::sqrt(5.0f)) * 0.5f; // z³ota proporcja
        std::vector<glm::vec3> verts = {
            glm::normalize(glm::vec3(-1,  t,  0)),
            glm::normalize(glm::vec3(1,  t,  0)),
            glm::normalize(glm::vec3(-1, -t,  0)),
            glm::normalize(glm::vec3(1, -t,  0)),

            glm::normalize(glm::vec3(0, -1,  t)),
            glm::normalize(glm::vec3(0,  1,  t)),
            glm::normalize(glm::vec3(0, -1, -t)),
            glm::normalize(glm::vec3(0,  1, -t)),

            glm::normalize(glm::vec3(t,  0, -1)),
            glm::normalize(glm::vec3(t,  0,  1)),
            glm::normalize(glm::vec3(-t,  0, -1)),
            glm::normalize(glm::vec3(-t,  0,  1)),
        };
        for (auto& v : verts) v *= radius; // na promieñ 'radius'

        // --- 20 œcian (CCW patrz¹c z zewn¹trz) ---
        std::vector<glm::ivec3> faces = {
            {0,11,5}, {0,5,1}, {0,1,7}, {0,7,10}, {0,10,11},
            {1,5,9},  {5,11,4}, {11,10,2}, {10,7,6}, {7,1,8},
            {3,9,4},  {3,4,2},  {3,2,6},  {3,6,8},  {3,8,9},
            {4,9,5},  {2,4,11}, {6,2,10}, {8,6,7},  {9,8,1}
        };

        // --- Subdivision ---
        for (int s = 0; s < SUBDIV; ++s) {
            std::unordered_map<uint64_t, int> cache;
            std::vector<glm::ivec3> new_faces;
            new_faces.reserve(faces.size() * 4);

            for (const auto& f : faces) {
                int i0 = f.x, i1 = f.y, i2 = f.z;

                int a = midpoint_index(i0, i1, verts, cache, radius);
                int b = midpoint_index(i1, i2, verts, cache, radius);
                int c = midpoint_index(i2, i0, verts, cache, radius);

                // 4 nowe trójk¹ty (CCW)
                new_faces.push_back({ i0, a, c });
                new_faces.push_back({ i1, b, a });
                new_faces.push_back({ i2, c, b });
                new_faces.push_back({ a,  b, c });
            }
            faces.swap(new_faces);
        }

        nodes.reserve(verts.size());
        glm::vec3 buf2 = glm::vec3(Source.src_x, Source.src_y, Source.src_z);
        glm::vec3 gSourceVel = Source.velocity;
        //const float audioE = getAtTime((gWinIdx)*window_ms / 1000.0f);
        const auto& wp = gWinPackets[gWinIdx];       // aktualne 5 ms
        int i = 0;
        for (const auto& p : verts) {
            node nd;
            nd.position = p;
            nd.nEmit = glm::normalize(p);
            nd.velocity = glm::normalize(p) * SOUND_V; // sta³a prêdkoœæ fali w oœrodku

            nd.position += buf2;
            nd.srcVel = gSourceVel;
            //nd.energy = audioE;
            //nd.tEmit = (gWinIdx)*gAudio.window_ms / 1000.0f; // zapisz czas emisji (sim-time)
            //nd.tEmit = (gWinIdx)*window_ms / 1000.0f;

            nd.energy = wp.amplitude;                 // u¿yj amplitudy okna
            nd.tEmit = wp.tEmit;                     // czas emisji okna (sek)
            //nd.winId = (int)gWinIdx;
            nd.seedId = i++;
            nodes.push_back(nd);
        }

        triangles.reserve(faces.size());
        for (const auto& f : faces) {
            triangles.push_back({ { f.x, f.y, f.z } });
            //break;
        }
        //buildSphereBuffers(/*dynamic=*/true);
        buildBuffersFor(nodes, triangles, gWaveGL, /*dynamic=*/true);

        // ...po zbudowaniu icosfery:
        gAvgEdgeLen = computeAvgEdgeLen(verts, faces);
        gBlindRadius = gBlindRings * gAvgEdgeLen;
        //gBlindRadius = 2.5f;

        first = false;
    }

    static int frameCount = 0;
    if ((frameCount % 8) == 0) {
        size_t budget = std::min<size_t>(4000, std::max<size_t>(1, triangles.size() / 10));
        budget = triangles.size(); //TO DO: MOZNA ZMIENIC NA WIEKSZE (W SENSIE ZMIENIC np. 4 -> 2)
        int threads = std::max(1u, std::thread::hardware_concurrency());

        if (refineIcosahedron_chunked_mt(0.1f, budget, threads)) {
            gMeshDirty = true; // przebuduj bufory tylko gdy zasz³a zmiana
        }
    }
    frameCount++;
    if (simulation_running)
    {
        updatePhysics(dt, Cube, Obstacle);
        //-----USUWANIE DOTKNIETYCH NODES-----
        pruneSlowNodes(/*minSpeed=*/getAtTime((gWinIdx)*window_ms / 1000.0f)); // m/s (dobierz)
    }

    // 1) odbuduj bufory tylko gdy trzeba
    if (gMeshDirty) {
        //buildSphereBuffers(/*dynamic=*/true);
        buildBuffersFor(nodes, triangles, gWaveGL, /*dynamic=*/true);
        gMeshDirty = false;
    }
    else {
        // 2) w przeciwnym razie tylko podmieñ pozycje
        //updateSpherePositions();
        updatePositionsFor(nodes, gWaveGL);
        updatePositionsFor(mic_nodes, gMicGL);
        updatePositionsFor(src_nodes, gSrcGL);
    }

    // 3) rysuj TYLKO z VBO/IBO
    //drawSphereWithBuffers();


    //kulka
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);
    //drawSphereWithBuffers();
    // FALA (bez offsetu — pozycje absolutne)
    drawMesh(gWaveGL);

    glEnable(GL_DEPTH_TEST);
    glm::vec3 micOffset = glm::vec3(Mic.mic_x, Mic.mic_y, Mic.mic_z);
    drawMesh(gMicGL, micOffset, /*fill*/{ 1.0f, 0.4f, 0.8f, 0.5f });
    glm::vec3 srcOffset = glm::vec3(Source.src_x, Source.src_y, Source.src_z);
    drawMesh(gSrcGL, srcOffset, /*fill*/{ 1.0f, 1.0f, 1.0f, 0.5f });
    //drawMicrophone();
    //drawSource();
    //PRZESZKODA

    glColor4f(0.2f, 0.5f, 0.2f, 0.7f);
    drawCuboidTransparentSorted(Obstacle);
    //glEnable(GL_DEPTH_TEST);
    // Basen
    glColor4f(0.0f, 0.0f, 1.0f, 0.1f);
    drawCuboidTransparentSorted(Cube);
}

int printOversizedTriangles(float maxArea) {
    // te same sta³e co w refine
    const float cuboidHalfWidth = Cube.width * 0.5f;
    const float cuboidHalfHeight = Cube.height * 0.5f;
    const float cuboidHalfDepth = Cube.depth * 0.5f;
    const float safetyMargin = 0.2f;

    // porównujemy |cross|^2 > (2*maxArea)^2
    const float thr2 = 4.0f * maxArea * maxArea;

    int count = 0;
    for (const auto& t : triangles) {
        int a = t.indices[0], b = t.indices[1], c = t.indices[2];
        const glm::vec3& pa = nodes[a].position;
        const glm::vec3& pb = nodes[b].position;
        const glm::vec3& pc = nodes[c].position;

        // pomijamy trójk¹ty blisko œciany (jak w refine)
        const bool nearWall =
            (std::fabs(pa.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pa.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pa.z) > cuboidHalfDepth - safetyMargin) ||
            (std::fabs(pb.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pb.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pb.z) > cuboidHalfDepth - safetyMargin) ||
            (std::fabs(pc.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pc.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pc.z) > cuboidHalfDepth - safetyMargin);

        if (nearWall) continue;

        // pole bez sqrt: |(pb-pa) x (pc-pa)|^2
        const glm::vec3 cr = glm::cross(pb - pa, pc - pa);
        const float cr2 = glm::dot(cr, cr);

        if (cr2 > thr2) ++count;
    }
    float procent = (float)count / (float)triangles.size() * 100.0f;
    std::cout << "Oversized (splittable) triangles: "
        << triangles.size() << "\n";
    return count;
}

static inline bool touchesMicrophone(const glm::vec3& p) {
    const glm::vec3 c(Mic.mic_x, Mic.mic_y, Mic.mic_z);
    const glm::vec3 d = p - c;
    // test sfery: odleg³oœæ^2 <= mic_radius^2
    return glm::dot(d, d) <= (mic_radius * mic_radius); // TO DO: trzeba sprawdzic dobre odleglosci, zeby mikrofon nie wykrywal za duzych kawalkow
}

// Usuwa wêz³y o |velocity| < minSpeed i remapuje trójk¹ty.
// Zwraca ile wêz³ów usuniêto.
int pruneSlowNodes(float minEnergy)
{
    if (nodes.empty()) return 0;

    const float thr2 = abs(minEnergy) / 10.0f;
    const size_t N = nodes.size();

    std::vector<int>  remap(N, -1);
    std::vector<node> kept;
    kept.reserve(N);

    bool only_one_read = false;
    //if (time_passed > 0.00095) doKill = true;
    int test = 0;
    for (size_t i = 0; i < N; ++i) {
        //nodes[i].energy = nodes[i].energy / nodes[i].path;
        const float E = abs(nodes[i].energy);
        //const glm::vec3 v = nodes[i].velocity;
        //const float v2 = glm::dot(v, v);

        bool EnergyEnough = (E > thr2);
        //if (nodes[i].bounces > 1) EnergyEnough = false;
        //bool EnergyEnough = true; // chwilowo do testow TO DO: ODKOMENTOWAC
        const bool hitMic = touchesMicrophone(nodes[i].position);

        if (hitMic) {
            const float win_s = window_ms / 1000.0f;         // 0.005 s
            float T = nodes[i].tEmit + time_passed; // czas przyjœcia pocz¹tku okna
            const int   wid = gWinIdx;              // ID okna z NODA (nie gWinIdx!)
            //if (nodes[i].bounces != 1) continue;

            // Je¿eli ten node jest czasowo st³umiony — nie loguj do CSV (ale „zu¿yj” go)
            if (T <= nodes[i].suppressUntilT) {
                only_one_read = true;
                continue;
            }
            //T = std::ceil((double)T * 100000.0) / 100000.0;

            /*if (nodes[i].bounces == 0 && T < prev_T2_time)
            {
                T = prev_T2_time;
            }*/

            ktore_odbicie++;
            // --- DOPPLER (rozci¹gniêcie/œciœniêcie osi czasu okna) ---
            //float predkosc = glm::length(nodes[i].velocity);
            const glm::vec3 n_hat = glm::normalize(nodes[i].velocity);
            const float     v_mic_proj = glm::dot(Mic.mic_velocity, n_hat) + 0.0415 * glm::length(Mic.mic_velocity);
            const float     v_mic_proj2 = glm::dot(n_hat, Mic.mic_velocity);
            const float v_src_proj = -glm::dot(nodes[i].srcVel, nodes[i].nEmit); // + gdy Ÿród³o jedzie w stronê propagacji

            //T2 = T + win_s * (1.0f + v_mic_proj / SOUND_V);
            float new_vel = nodes[i].path / time_passed;
            //T2 = T + win_s * ((new_vel + v_src_proj) / (new_vel - v_mic_proj));
            float proba = T + win_s * ((SOUND_V + v_src_proj) / (SOUND_V - v_mic_proj));
            //T2 = T + win_s * ((SOUND_V + v_src_proj) / (SOUND_V - v_mic_proj)) - 10*dt*(nodes[i].bounces+1);
            //T2 = T + win_s * ((new_vel + v_src_proj) / (new_vel - v_mic_proj));
            float vel_test = glm::length(nodes[i].velocity);
            //float vel_test2 = nodes[i].velocity.length();
            T2 = T + win_s * nodes[i].doppler * ((vel_test + v_src_proj) / (vel_test - v_mic_proj));
            //T2 = std::ceil((double)T2 * 100000.0) / 100000.0;

            //if (!rewind_punkt)
            //{
            //    glm::vec3 droga = glm::vec3(Mic.mic_x, Mic.mic_y, Mic.mic_z) - Mic.rewind_point;
            //    float droga2 = glm::length(droga);
            //    //T2 = T + (win_s * (predkosc + v_src_proj) - droga2) / (predkosc - v_mic_proj);
            //    //T2 = T + nodes[i].doppler * (win_s * (SOUND_V + v_src_proj) - droga2) / (SOUND_V - v_mic_proj);
            //}
            //T2 = (float)(std::trunc((double)T2 * 100000.0) / 100000.0);
            //prev_T2_time = T2;

            //float doppler = ((SOUND_V + v_src_proj) / (SOUND_V - v_mic_proj));
            //T2 = T + win_s * 1 / doppler;

            // --- CLAMP czasu startu okna, gdy "oddalamy siê" (scale < 1) ---
            if (T < prev_T2_time) {
                //T = prev_T2_time;               // wymuœ ci¹g³oœæ z poprzednim T2
                //T2 = T + win_s * doppler;          // zachowaj ten sam wspó³czynnik Dopplera
                int siema = 3;
            }



            float           scale = (T2 - T) / win_s;
            if (scale <= 0.0f) scale = 1e-4f;

            // Logowanie pe³nego okna 5 ms (czas + wartoœci), amplituda dopasowana do noda
            if (wid >= 0 && wid < (int)gWinPackets.size()) {
                const auto& wp = gWinPackets[wid];
                const float  Awin = (wp.amplitude != 0.0f) ? wp.amplitude : 1e-12f;
                //const float  scaleAmp = nodes[i].energy / Awin;
                // --- przed pêtl¹ po j ---
                constexpr float TS = 1e-5f;   // odstêp czasu jaki chcemy wymuszaæ
                constexpr float EPS = 2e-7f;   // ma³a tolerancja na b³êdy binarne
                float t_last = std::numeric_limits<float>::quiet_NaN();
                float v_last = 0.0f;

                // --- korekta amplitudy wewn¹trz okna 5 ms (zmiana zasiêgu r(t)) ---
                const float win_s = window_ms / 1000.0f;      // 0.005 s
                const float r0 = std::max(nodes[i].path, R0_ATTEN);


                // U Ciebie: v_mic_proj = dot(Mic.mic_velocity, n_hat)
                //           v_src_proj = -dot(nodes[i].srcVel, nodes[i].nEmit)  // ju¿ policzone wy¿ej
                const float rdot = v_mic_proj + v_src_proj;


                const float PATH_SCALE = 0.00001f;


                for (size_t j = 0; j < wp.times.size(); ++j) {
                    float t_abs = T + scale * wp.times[j];     // [T, T2] po Dopplerz


                    const float  scaleAmp = nodes[i].energy / Awin;

                    // czas wzglêdny w *odbiorze*, wzglêdem œrodka okna (¿eby nie biasowaæ ca³oœci)
                    float t_loc_rec = scale * (wp.times[j] - 0.5f * win_s);

                    // efektywny zasiêg dla tej próbki (spójny z Twoim units: PATH_SCALE)
                    float r_eff = std::max(r0 + PATH_SCALE * rdot * t_loc_rec, R0_ATTEN);


                    float intra = r0 / r_eff;

                    // dotychczas: float val = wp.values[j] * scaleAmp;
                    float val = (wp.values[j] * scaleAmp) * intra;

                    //float val = wp.values[j] * scaleAmp;     // amplituda po odbiciach
                    //t_abs += nodes[i].bounces * 0.2f;          // jeœli chcesz zachowaæ offset odbiæ
                    //t_abs += gWinIdx * 0.1f;
                    //t_abs += ktore_odbicie * 0.05f;
                    if (t_abs >= 0.00977f)
                    {
                        t_abs = t_abs;
                        int i = 2;
                    }
                    float t_out = q5(t_abs);
                    if (t_out == t_last) continue;
                    // jeœli mamy poprzedni punkt i luka jest wiêksza ni¿ TS, to do³ó¿ brakuj¹ce kroki
                    if (!std::isnan(t_last)) {
                        float gap = t_out - t_last;

                        if (gap > TS + EPS) {
                            // ile kroków TS mieœci siê w luce? (ile „szczelin” + 1)
                            int steps = (int)std::floor((gap + EPS) / TS);

                            // wstaw brakuj¹ce punkty: t_last + s*TS, s = 1 .. steps-1
                            for (int s = 1; s < steps; ++s) {
                                float tf = t_last + s * TS;
                                // liniowa interpolacja wartoœci miêdzy (t_last,v_last) a (t_out,val)
                                float a = (tf - t_last) / gap;
                                float vf = v_last + a * (val - v_last);
                                logMicHit(q5(tf), vf);
                            }
                        }
                    }

                    // zapisz punkt bie¿¹cy
                    logMicHit(q5(t_out), val);

                    // aktualizuj „ostatni zapisany”
                    t_last = t_out;
                    v_last = val;
                }
            }
            else {
                // Fallback: gdyby brakowa³o pakietu okna
                logMicHit(T, nodes[i].energy);
            }

            // === OŒLEPIANIE WY£¥CZNIE TEGO SAMEGO CZO£A (lokalny klaster) ===
            // (1) horyzont czasowy t³umienia
            // Wersja sta³a:
            const float blind_until = T + 0.0015f;
            // Wersja „fizyczna” (zale¿na od geometrii i prêdkoœci wzglêdnej) — jeœli wolisz:
            // float c_rel = std::max(SOUND_V - v_mic_proj, 0.1f * SOUND_V);
            // const float blind_until = T + (gBlindRadius / c_rel);

            // (2) parametry filtrów klastra
            const glm::vec3 center = nodes[i].position;                // pozycja trafionego noda
            const glm::vec3 micC = glm::vec3(Mic.mic_x, Mic.mic_y, Mic.mic_z);   // <-- PODMIEÑ NA SWOJE POLA!
            const float     cosMax = std::cos(glm::radians(20.0f));    // max ró¿nica kierunku (~20°)
            const float     s_i = glm::dot(micC - nodes[i].position, n_hat); // proxy czasu dojœcia

            for (size_t j = 0; j < nodes.size(); ++j) {
                if (j == i) continue;

                // To samo okno emisji (to samo czo³o czasowo)
                if (gWinIdx != wid) continue;

                // Zgodny kierunek propagacji (wyklucz inne odbicia)
                const glm::vec3 nj_hat = glm::normalize(nodes[j].velocity);
                if (glm::dot(n_hat, nj_hat) < cosMax) continue;

                // Node musi zbli¿aæ siê do mikrofonu
                const float s_j = glm::dot(micC - nodes[j].position, nj_hat);
                if (s_j <= 0.0f) continue;

                // Zbli¿ony przewidywany czas dojœcia (ró¿nica projekcji)
                if (std::fabs(s_j - s_i) > 1.5f * gBlindRadius) continue;

                // Blisko w przestrzeni (ok. „3 trójk¹ty”)
                if (glm::length(nodes[j].position - center) > gBlindRadius) continue;

                // Ustaw t³umienie do blind_until — zarówno w nodes[], jak i ewentualnie w kept[]
                if (remap[j] != -1) {
                    const int kj = remap[j];
                    kept[kj].suppressUntilT = std::max(kept[kj].suppressUntilT, blind_until);
                }
                else {
                    nodes[j].suppressUntilT = std::max(nodes[j].suppressUntilT, blind_until);
                }
            }

            doKill = true;
            only_one_read = true;   // jeœli u Ciebie „zu¿ywa” noda po trafieniu
            // (opcjonalnie) debug:
            std::cout << "win=" << wid << "  T=" << T << "\n";
            //std::cout << ktore_odbicie << std::endl;
            continue;
        }


        if (EnergyEnough) {
            remap[i] = (int)kept.size();
            node cp = nodes[i];
            // NA WSZELKI WYPADEK zsynchronizuj suppress z ewentualnymi pendingami:
            // (jeœli wczeœniej ktoœ dopisa³ do nodes[i], to ju¿ to masz; ten krok tylko wzmacnia spójnoœæ)
            // cp.suppressUntilT = std::max(cp.suppressUntilT, ???); // niepotrzebne przy remap-aktualizacji wy¿ej
            kept.push_back(cp);
        }
    }
    // pamiêtaj najnowszy koniec bie¿¹cej ramki
    prev_T2_time = T2;

    // Triangles: zostaw tylko te, które przetrwa³y
    std::vector<Triangle> newTris;
    newTris.reserve(triangles.size());
    for (const auto& t : triangles) {
        int na = remap[t.indices[0]];
        int nb = remap[t.indices[1]];
        int nc = remap[t.indices[2]];
        if (na >= 0 && nb >= 0 && nc >= 0) {
            newTris.push_back({ { na, nb, nc } });
        }
    }

    nodes.swap(kept);
    triangles.swap(newTris);

    // Uniewa¿nij cache/rafineriê i poproœ o rebuild VBO
    gRefineCursor = 0;
    gEdgeMidCache.clear();
    gMeshDirty = true;

    if (doKill) killAllNodes(); // do testow
    // zgasniecie fali
    const bool windowDied = nodes.empty();

    if (windowDied) {
        if (!beginNextWindow() or gWinIdx == windows_number) {
            writeMicCsv(write_file);
            resetMicEvents();
            simulation_running = false;
            std::cout << "KONIEC SYMULACJI" << std::endl;
        }
    }

    return (int)(N - nodes.size());
}

void drawCuboidTransparentSorted(struct Cuboid_dimensions temp_Cube) {
    float halfWidth = temp_Cube.width / 2.0f;
    float halfHeight = temp_Cube.height / 2.0f;
    float halfDepth = temp_Cube.depth / 2.0f;

    glm::vec3 vertices[] = {
         { -halfWidth + temp_Cube.x_offset, -halfHeight + temp_Cube.y_offset,  halfDepth + temp_Cube.z_offset},
         {  halfWidth + temp_Cube.x_offset, -halfHeight + temp_Cube.y_offset,  halfDepth + temp_Cube.z_offset},
         {  halfWidth + temp_Cube.x_offset,  halfHeight + temp_Cube.y_offset,  halfDepth + temp_Cube.z_offset},
         { -halfWidth + temp_Cube.x_offset,  halfHeight + temp_Cube.y_offset,  halfDepth + temp_Cube.z_offset},
         { -halfWidth + temp_Cube.x_offset, -halfHeight + temp_Cube.y_offset, -halfDepth + temp_Cube.z_offset},
         {  halfWidth + temp_Cube.x_offset, -halfHeight + temp_Cube.y_offset, -halfDepth + temp_Cube.z_offset},
         {  halfWidth + temp_Cube.x_offset,  halfHeight + temp_Cube.y_offset, -halfDepth + temp_Cube.z_offset},
         { -halfWidth + temp_Cube.x_offset,  halfHeight + temp_Cube.y_offset, -halfDepth + temp_Cube.z_offset}
    };

    struct Face {
        int indices[4];
        float distance;
    };

    Face faces[6] = {
        {{0, 1, 2, 3}, 0}, // przód
        {{4, 5, 6, 7}, 0}, // ty³
        {{0, 3, 7, 4}, 0}, // lewo
        {{1, 2, 6, 5}, 0}, // prawo
        {{0, 1, 5, 4}, 0}, // dó³
        {{2, 3, 7, 6}, 0}  // góra
    };

    // Oblicz odleg³oœci œcian od kamery
    for (int i = 0; i < 6; ++i) {
        glm::vec3 center(0.0f);
        for (int j = 0; j < 4; ++j) {
            center += vertices[faces[i].indices[j]];
        }
        center /= 4.0f;
        faces[i].distance = glm::length(center - cameraPos);
    }

    // Sortuj œciany od najdalszej do najbli¿szej
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            if (faces[i].distance < faces[j].distance) {
                std::swap(faces[i], faces[j]);
            }
        }
    }

    // Rysuj œciany w odpowiedniej kolejnoœci
    glBegin(GL_QUADS);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            glVertex3f(vertices[faces[i].indices[j]].x,
                vertices[faces[i].indices[j]].y,
                vertices[faces[i].indices[j]].z);
        }
    }
    glEnd();
}

// Zabija wszystkie nody (i trójk¹ty) – do debugu.
void killAllNodes()
{
    nodes.clear();
    triangles.clear();

    // wyczyœæ stan rafinowania/mesh
    gRefineCursor = 0;
    gEdgeMidCache.clear();
    gMeshDirty = true;

    //gIndexCount = 0;

    // opcjonalnie: od razu „opró¿nij” bufory GPU
    //if (gVBO) { glBindBuffer(GL_ARRAY_BUFFER, gVBO); glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW); glBindBuffer(GL_ARRAY_BUFFER, 0); }
    //if (gIBO) { glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIBO); glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW); glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); }
}

void buildBuffersFor(const std::vector<node>& verts,
    const std::vector<Triangle>& tris,
    MeshGL& m,
    bool dynamic)
{
    if (!GLAD_GL_VERSION_1_5) { std::cerr << "VBO niedostêpne (GL < 1.5)\n"; return; }
    m.dynamic = dynamic;

    if (!m.vbo) glGenBuffers(1, &m.vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m.vbo);

    const GLsizeiptr vbSize = (GLsizeiptr)(verts.size() * sizeof(node));
    const GLenum usage = dynamic ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW;
    glBufferData(GL_ARRAY_BUFFER, vbSize, verts.empty() ? nullptr : (const void*)verts.data(), usage);

    if (!m.ibo) glGenBuffers(1, &m.ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ibo);

    const GLsizeiptr ibSize = (GLsizeiptr)(tris.size() * 3 * sizeof(unsigned int));
    const void* idxSrc = tris.empty() ? nullptr : (const void*)&tris[0].indices[0];
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, ibSize, idxSrc, GL_STATIC_DRAW);

    m.indexCount = (GLsizei)(tris.size() * 3);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void updatePositionsFor(const std::vector<node>& verts, MeshGL& m)
{
    if (!m.vbo) return;
    glBindBuffer(GL_ARRAY_BUFFER, m.vbo);

    const GLsizeiptr vbSize = (GLsizeiptr)(verts.size() * sizeof(node));
    glBufferData(GL_ARRAY_BUFFER, vbSize, nullptr, GL_DYNAMIC_DRAW);            // orphan
    if (!verts.empty())
        glBufferSubData(GL_ARRAY_BUFFER, 0, vbSize, (const void*)verts.data()); // upload

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void drawMesh(const MeshGL& m,
    const glm::vec3& offset,
    const glm::vec4& fill,
    const glm::vec4& wire,
    float wireWidth)
{
    if (!m.vbo || !m.ibo || m.indexCount == 0) return;

    glPushMatrix();
    if (offset.x != 0 || offset.y != 0 || offset.z != 0)
        glTranslatef(offset.x, offset.y, offset.z);

    glBindBuffer(GL_ARRAY_BUFFER, m.vbo);
    glEnableClientState(GL_VERTEX_ARRAY); // compat profile
    glVertexPointer(3, GL_FLOAT, (GLsizei)sizeof(node), (const void*)offsetof(node, position));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ibo);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor4f(fill.r, fill.g, fill.b, fill.a);
    glDrawElements(GL_TRIANGLES, m.indexCount, GL_UNSIGNED_INT, (void*)0);

    glDisable(GL_POLYGON_OFFSET_FILL);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth(wireWidth);
    glColor4f(wire.r, wire.g, wire.b, wire.a);
    glDrawElements(GL_TRIANGLES, m.indexCount, GL_UNSIGNED_INT, (void*)0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glPopMatrix();
}

// Zwraca sferê wokó³ (0,0,0) o promieniu 'radius' i gêstoœci 'subdiv'
static void buildIcosphereNodesTris(float radius, int subdiv,
    std::vector<node>& out_nodes,
    std::vector<Triangle>& out_tris)
{
    // --- 12 wierzcho³ków icosa ---
    const float t = (1.0f + std::sqrt(5.0f)) * 0.5f;
    std::vector<glm::vec3> verts = {
        glm::normalize(glm::vec3(-1,  t,  0)),
        glm::normalize(glm::vec3(1,  t,  0)),
        glm::normalize(glm::vec3(-1, -t,  0)),
        glm::normalize(glm::vec3(1, -t,  0)),

        glm::normalize(glm::vec3(0, -1,  t)),
        glm::normalize(glm::vec3(0,  1,  t)),
        glm::normalize(glm::vec3(0, -1, -t)),
        glm::normalize(glm::vec3(0,  1, -t)),

        glm::normalize(glm::vec3(t,  0, -1)),
        glm::normalize(glm::vec3(t,  0,  1)),
        glm::normalize(glm::vec3(-t,  0, -1)),
        glm::normalize(glm::vec3(-t,  0,  1)),
    };
    for (auto& v : verts) v *= radius;

    std::vector<glm::ivec3> faces = {
        {0,11,5}, {0,5,1}, {0,1,7}, {0,7,10}, {0,10,11},
        {1,5,9},  {5,11,4}, {11,10,2}, {10,7,6}, {7,1,8},
        {3,9,4},  {3,4,2},  {3,2,6},  {3,6,8},  {3,8,9},
        {4,9,5},  {2,4,11}, {6,2,10}, {8,6,7},  {9,8,1}
    };

    // --- Subdivision z cache krawêdzi ---
    for (int s = 0; s < subdiv; ++s) {
        std::unordered_map<uint64_t, int> cache;
        std::vector<glm::ivec3> new_faces;
        new_faces.reserve(faces.size() * 4);

        for (const auto& f : faces) {
            int i0 = f.x, i1 = f.y, i2 = f.z;
            int a = midpoint_index(i0, i1, verts, cache, radius);
            int b = midpoint_index(i1, i2, verts, cache, radius);
            int c = midpoint_index(i2, i0, verts, cache, radius);

            new_faces.push_back({ i0, a, c });
            new_faces.push_back({ i1, b, a });
            new_faces.push_back({ i2, c, b });
            new_faces.push_back({ a,  b, c });
        }
        faces.swap(new_faces);
    }

    // --- Konwersja do Twoich typów ---
    out_nodes.clear();
    out_nodes.reserve(verts.size());
    for (const auto& p : verts) {
        node nd;
        nd.position = p;                 // MODEL-SPACE (bez translacji!)
        nd.velocity = glm::vec3(0);      // statyczne
        nd.energy = 0.0f;
        out_nodes.push_back(nd);
    }

    out_tris.clear();
    out_tris.reserve(faces.size());
    for (const auto& f : faces) {
        out_tris.push_back({ { f.x, f.y, f.z } });
    }
}

