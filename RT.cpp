#include <iostream>   // cout, cerr
#include <fstream>    // ofstream
#include <vector>     // vector
#include <cmath>      // sqrt, pow, sin, cos, M_PI
#include <limits>     // numeric_limits
#include <iomanip>    // setw, setfill
#include <sstream>    // ostringstream

using namespace std;

const float INF = numeric_limits<float>::infinity();
const int Cw = 400, Ch = 400;
const float Vw = 1, Vh = 1;
const float d = 1;
const int TOTAL_FRAMES = 60;


class Color {
public:
    int r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(int r, int g, int b) : r(r), g(g), b(b) {}
};

const Color BACKGROUND_COLOR(0, 0, 0);


class Vector3D {
public:
    float x, y, z;
    
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(float x, float y, float z) : x(x), y(y), z(z) {}
    
    float length() const {
        return sqrt(x*x + y*y + z*z);
    }
    
    Vector3D normalize() const {
        float len = length();
        return Vector3D(x/len, y/len, z/len);
    }
    
    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }
    
    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }
    
    Vector3D operator*(float scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }
    
    Vector3D operator/(float scalar) const {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }
};

float dot(const Vector3D& a, const Vector3D& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector3D operator*(float scalar, const Vector3D& v) {
    return v * scalar;
}

// Refleja el vector R sobre la normal N
Vector3D reflectRay(const Vector3D& R, const Vector3D& N) {
    return 2 * dot(N, R) * N - R;
}


class Sphere {
public:
    Vector3D center;
    float radius;
    Color color;
    float specular;
    float reflection;
    
    Sphere(Vector3D center, float radius, Color color, float specular, float reflection)
        : center(center), radius(radius), color(color), specular(specular), reflection(reflection) {}
};


class Plane {
public:
    Vector3D normal;
    float distance;
    Color color;
    float specular;
    float reflection;
    
    Plane(Vector3D normal, float distance, Color color, float specular, float reflection)
        : normal(normal.normalize()), distance(distance), color(color), specular(specular), reflection(reflection) {}
};


class Light {
public:
    string type;
    float intensity;
    Vector3D position;
    Vector3D direction;
    
    Light(string type, float intensity, Vector3D position = Vector3D(), Vector3D direction = Vector3D())
        : type(type), intensity(intensity), position(position), direction(direction) {}
};

class Scene {
public:
    vector<Sphere> spheres;
    vector<Plane> planes;
    vector<Light> lights;
    
    Scene() {
        // Inicializa un plano
        planes.push_back(Plane(Vector3D(0, 1, 0), -1.0f, Color(255, 255, 255), 500, 0.2f));
        // Configurar las luces
        lights.push_back(Light("ambient", 0.2f));
        lights.push_back(Light("point", 0.6f, Vector3D(2, 1, 0)));
        lights.push_back(Light("directional", 0.2f, Vector3D(), Vector3D(1, 4, 4)));
    }
    
    // Para actualizar la escena para la animacion recalculando las posiciones de las esferas
    void update(int frame) {
        spheres.clear();
        float t = frame / float(TOTAL_FRAMES);
        float height1 = 0.5f + 0.5f * sin(t * 2 * M_PI);
        float height2 = -0.5f - 0.5f * sin(t * 2 * M_PI);
        Vector3D pos1(0.7f, height1, 3.0f);
        Vector3D pos2(-0.7f, height2, 3.0f);
        spheres.push_back(Sphere(pos1, 0.5f, Color(236, 124, 38), 1000, 0));
        spheres.push_back(Sphere(pos2, 0.5f, Color(123, 31, 25), 500, 0.2f));
    }
};

Scene scene;

pair<float, float> intersectRaySphere(const Vector3D& O, const Vector3D& D, const Sphere& sphere) {
    Vector3D CO = O - sphere.center;
    float a = dot(D, D);
    if (a == 0) return make_pair(INF, INF);
    
    float b = 2 * dot(CO, D);
    float c = dot(CO, CO) - sphere.radius * sphere.radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return make_pair(INF, INF);
    
    float sqrtDisc = sqrt(discriminant);
    float t1 = (-b + sqrtDisc) / (2 * a);
    float t2 = (-b - sqrtDisc) / (2 * a);
    
    return make_pair(t1, t2);
}

float intersectRayPlane(const Vector3D& O, const Vector3D& D, const Plane& plane) {
    float denom = dot(plane.normal, D);
    if (fabs(denom) < 1e-6)
        return INF;
    
    float t = (plane.distance - dot(plane.normal, O)) / denom;
    return (t >= 0) ? t : INF;
}

pair<Sphere*, float> closestSphereIntersection(const Vector3D& O, const Vector3D& D, float tMin, float tMax) {
    float closestT = INF;
    Sphere* closestSphere = nullptr;
    
    for (auto& sphere : scene.spheres) {
        auto [t1, t2] = intersectRaySphere(O, D, sphere);
        if (t1 > tMin && t1 < tMax && t1 < closestT) {
            closestT = t1;
            closestSphere = &sphere;
        }
        if (t2 > tMin && t2 < tMax && t2 < closestT) {
            closestT = t2;
            closestSphere = &sphere;
        }
    }
    
    return make_pair(closestSphere, closestT);
}

pair<Plane*, float> closestPlaneIntersection(const Vector3D& O, const Vector3D& D, float tMin, float tMax) {
    float closestT = INF;
    Plane* closestPlane = nullptr;
    
    for (auto& plane : scene.planes) {
        float t = intersectRayPlane(O, D, plane);
        if (t > tMin && t < tMax && t < closestT) {
            closestT = t;
            closestPlane = &plane;
        }
    }
    
    return make_pair(closestPlane, closestT);
}

// Devuelve la interseccion mas cercana, ya sea una esfera o un plano
pair<void*, float> closestIntersection(const Vector3D& O, const Vector3D& D, float tMin, float tMax) {
    auto [spherePtr, sphereT] = closestSphereIntersection(O, D, tMin, tMax);
    auto [planePtr, planeT] = closestPlaneIntersection(O, D, tMin, tMax);
    
    if (sphereT < planeT)
        return make_pair(static_cast<void*>(spherePtr), sphereT);
    else
        return make_pair(static_cast<void*>(planePtr), planeT);
}

float computeLighting(const Vector3D& P, const Vector3D& N, const Vector3D& V, float specular) {
    float intensity = 0.0f;
    for (const auto& light : scene.lights) {
        if (light.type == "ambient") {
            intensity += light.intensity;
        } else {
            Vector3D L = (light.type == "point") ? light.position - P : light.direction;
            L = L.normalize();
            // Check for shadows
            auto [shadowObj, shadowT] = closestIntersection(P, L, 0.001f, INF);
            if (shadowObj != nullptr)
                continue;
            float nDotL = dot(N, L);
            if (nDotL > 0)
                intensity += light.intensity * nDotL / (N.length() * L.length());
            if (specular > 0) {
                Vector3D R = reflectRay(L, N);
                float rDotV = dot(R, V);
                if (rDotV > 0)
                    intensity += light.intensity * pow(rDotV, specular);
            }
        }
    }
    return intensity;
}

Color traceRay(const Vector3D& O, const Vector3D& D, float tMin, float tMax, int recursionDepth) {
    auto [closestObj, closestT] = closestIntersection(O, D, tMin, tMax);
    if (!closestObj)
        return BACKGROUND_COLOR;
    
    Vector3D P = O + D * closestT;
    Vector3D N;
    Color color;
    float specular, reflection;
    
    if (Sphere* sphere = dynamic_cast<Sphere*>(static_cast<Sphere*>(closestObj))) {
        N = (P - sphere->center).normalize();
        color = sphere->color;
        specular = sphere->specular;
        reflection = sphere->reflection;
    } else if (Plane* plane = dynamic_cast<Plane*>(static_cast<Plane*>(closestObj))) {
        N = plane->normal;
        color = plane->color;
        specular = plane->specular;
        reflection = plane->reflection;
    } else {
        return BACKGROUND_COLOR;
    }
    
    Vector3D V = (-1.0f * D).normalize();
    float lightIntensity = computeLighting(P, N, V, specular);
    Color localColor(
        min(255, static_cast<int>(color.r * lightIntensity)),
        min(255, static_cast<int>(color.g * lightIntensity)),
        min(255, static_cast<int>(color.b * lightIntensity))
    );
    
    if (recursionDepth <= 0 || reflection <= 0)
        return localColor;
    
    Vector3D R = reflectRay(V, N);
    Color reflectedColor = traceRay(P, R, 0.001f, INF, recursionDepth - 1);
    
    return Color(
        static_cast<int>((1 - reflection) * localColor.r + reflection * reflectedColor.r),
        static_cast<int>((1 - reflection) * localColor.g + reflection * reflectedColor.g),
        static_cast<int>((1 - reflection) * localColor.b + reflection * reflectedColor.b)
    );
}

// Convierte las coordenadas del canvas en coordenadas del viewport
Vector3D canvasToViewport(float x, float y) {
    return Vector3D(x * Vw / Cw, y * Vh / Ch, d);
}

void writePPM(const vector<Color>& pixels, const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    file << "P3\n" << Cw << " " << Ch << "\n255\n";
    for (const auto& pixel : pixels)
        file << pixel.r << " " << pixel.g << " " << pixel.b << "\n";
}

void renderFrame(int frame) {
    vector<Color> pixels(Cw * Ch);
    Vector3D O(0, 0, 0);
    scene.update(frame);
    for (int y = 0; y < Ch; y++) {
        for (int x = 0; x < Cw; x++) {
            float canvasX = x - Cw/2;
            float canvasY = Ch/2 - y;
            Vector3D D = canvasToViewport(canvasX, canvasY);
            pixels[y * Cw + x] = traceRay(O, D, 1, INF, 3);
        }
    }
    ostringstream filename;
    filename << "frames/frame_" << setw(4) << setfill('0') << frame << ".ppm";
    writePPM(pixels, filename.str());
    cout << "Rendered frame " << frame << " of " << TOTAL_FRAMES << endl;
}

int main() {
    system("mkdir  frames");
    for (int frame = 1; frame < TOTAL_FRAMES; frame++) {
        renderFrame(frame);
    }
    cout << "Combining frames into a video." << endl;
    system("ffmpeg -y -framerate 30 -i frames/frame_%04d.ppm -c:v libx264 -pix_fmt yuv420p output.mp4");
    cout << "Rendering complete. Video saved to output.mp4" << endl;
    return 0;
}
