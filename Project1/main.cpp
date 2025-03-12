#include <GLFW/glfw3.h>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#define M_PI 3.1415

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void renderScene();
// do ogarniecia
void drawCuboid(float width, float height, float depth);
void drawOrangeSphereWithBlackEdges(float radius, int segments);
void drawIcosahedron(float radius);

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

struct Vec3 {
    float x, y, z;
    Vec3 operator+(const Vec3& other) const { return { x + other.x, y + other.y, z + other.z }; }
    Vec3 operator-(const Vec3& other) const { return { x - other.x, y - other.y, z - other.z }; }
    Vec3 operator*(float scalar) const { return { x * scalar, y * scalar, z * scalar }; }
    Vec3 cross(const Vec3& other) const { return { y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x }; }
    Vec3 normalize() const {
        float len = std::sqrt(x * x + y * y + z * z);
        return { x / len, y / len, z / len };
    }
};


Vec3 cameraPos = { 0.0f, 0.0f, 3.0f };
Vec3 cameraFront = { 0.0f, 0.0f, -1.0f };
Vec3 cameraUp = { 0.0f, 1.0f, 0.0f };

bool firstMouse = true;
bool mousePressed = false;
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
float yaw = -90.0f;
float pitch = 0.0f;
float fov = 45.0f;

float deltaTime = 0.0f;
float lastFrame = 0.0f;

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
    // Zmiana odleg³oœci kamery za pomoc¹ scrolla
    cameraPos = cameraPos + cameraFront * static_cast<float>(yoffset) * 0.5f;
}

int main()
{
    if (!glfwInit())
        return -1;

    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Camera Scene", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback); // Dodanie obs³ugi scrolla
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    //glEnable(GL_DEPTH_TEST); // wazne ???? 
    //glEnable(GL_DEPTH_TEST || GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);

        glClearColor(0.5f, 0.5f, 0.5f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 projection = glm::perspective(glm::radians(fov), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = glm::lookAt(glm::vec3(cameraPos.x, cameraPos.y, cameraPos.z),
            glm::vec3(cameraPos.x + cameraFront.x, cameraPos.y + cameraFront.y, cameraPos.z + cameraFront.z),
            glm::vec3(cameraUp.x, cameraUp.y, cameraUp.z));

        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(&projection[0][0]);

        glMatrixMode(GL_MODELVIEW);
        glLoadMatrixf(&view[0][0]);

        renderScene();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

void processInput(GLFWwindow* window)
{
    float cameraSpeed = 2.5f * deltaTime;

    // Zmiana wysokoœci kamery za pomoc¹ klawiszy W i S
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos.y += cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos.y -= cameraSpeed;

    // Poruszanie kamer¹ w p³aszczyŸnie XZ za pomoc¹ klawiszy A i D
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos = cameraPos - cameraFront.cross(cameraUp).normalize() * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos = cameraPos + cameraFront.cross(cameraUp).normalize() * cameraSpeed;
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

    Vec3 front;
    front.x = cos(yaw * M_PI / 180.0) * cos(pitch * M_PI / 180.0);
    front.y = sin(pitch * M_PI / 180.0);
    front.z = sin(yaw * M_PI / 180.0) * cos(pitch * M_PI / 180.0);
    cameraFront = front.normalize();
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void renderScene()
{
    glEnable(GL_BLEND); // W³¹cz blending (mieszanie kolorów)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Ustaw funkcjê mieszania
    //glEnable(GL_DEPTH_TEST); // moze wazne

    //glBegin(GL_QUADS);
    // 
    //glColor4f(0.0f, 1.0f, 0.0f, 0.3f); // Kolor niebieski (alfa = 1.0, nieprzezroczysty)
    //glVertex3f(-0.5f, -0.5f, 0.0f);
    //glVertex3f(0.5f, -0.5f, 0.0f);
    //glVertex3f(0.5f, 0.5f, 0.0f);
    //glVertex3f(-0.5f, 0.5f, 0.0f);
    //
    //glColor4f(0.0f, 0.0f, 1.0f, 0.3f); // Kolor niebieski (alfa = 1.0, nieprzezroczysty)
    //glVertex3f(-0.5f, -0.5f, 0.5f);
    //glVertex3f(0.5f, -0.5f, 0.5f);
    //glVertex3f(0.5f, 0.5f, 0.5f);
    //glVertex3f(-0.5f, 0.5f, 0.5f);

    //glColor4f(1.0f, 0.0f, 0.0f, 0.3f);
    //glVertex3f(-0.5f, -0.5f, -0.5f);
    //glVertex3f(0.5f, -0.5f, -0.5f);
    //glVertex3f(0.5f, 0.5f, -0.5f);
    //glVertex3f(-0.5f, 0.5f, -0.5f); 
    //glEnd();

    drawIcosahedron(0.5f);
    glColor4f(0.0f, 0.0f, 1.0f, 0.3f);
    drawCuboid(8.0f, 6.0f, 6.0f);

    //glDisable(GL_BLEND); // Wy³¹cz blending po zakoñczeniu rysowania
    //drawIcosahedron(0.5f);
    //glDisable(GL_DEPTH_TEST);
}
void drawCuboid(float width, float height, float depth) {
    float halfWidth = width / 2.0f;
    float halfHeight = height / 2.0f;
    float halfDepth = depth / 2.0f;

    // Definiowanie wierzcho³ków prostopad³oœcianu
    glm::vec3 vertices[] = {
        // Przednia œciana
        { -halfWidth, -halfHeight,  halfDepth }, // 0
        {  halfWidth, -halfHeight,  halfDepth }, // 1
        {  halfWidth,  halfHeight,  halfDepth }, // 2
        { -halfWidth,  halfHeight,  halfDepth }, // 3

        // Tylna œciana
        { -halfWidth, -halfHeight, -halfDepth }, // 4
        {  halfWidth, -halfHeight, -halfDepth }, // 5
        {  halfWidth,  halfHeight, -halfDepth }, // 6
        { -halfWidth,  halfHeight, -halfDepth }  // 7
    };

    // Indeksy wierzcho³ków dla œcian (czworok¹tów)
    int faces[6][4] = {
        {0, 1, 2, 3}, // Przednia œciana
        {4, 5, 6, 7}, // Tylna œciana
        {0, 3, 7, 4}, // Lewa œciana
        {1, 2, 6, 5}, // Prawa œciana
        {0, 1, 5, 4}, // Dolna œciana
        {2, 3, 7, 6}  // Górna œciana
    };

    // Rysowanie œcian (czworok¹tów)
    glBegin(GL_QUADS);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            glVertex3f(vertices[faces[i][j]].x, vertices[faces[i][j]].y, vertices[faces[i][j]].z);
        }
    }
    glEnd();

    // Rysowanie krawêdzi (linii)
    glColor3f(0.0f, 0.0f, 0.0f); // Ustaw kolor na czarny
    glBegin(GL_LINES);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            int current = faces[i][j];
            int next = faces[i][(j + 1) % 4]; // Nastêpny wierzcho³ek w œcianie
            glVertex3f(vertices[current].x, vertices[current].y, vertices[current].z);
            glVertex3f(vertices[next].x, vertices[next].y, vertices[next].z);
        }
    }
    glEnd();
}
void drawOrangeSphereWithBlackEdges(float radius, int segments) {
    const float PI = 3.14159265358979323846f;

    for (int i = 0; i < segments; ++i) {
        float theta1 = i * PI / segments;       // K¹t theta (pionowy)
        float theta2 = (i + 1) * PI / segments;

        for (int j = 0; j < segments; ++j) {
            float phi1 = j * 2 * PI / segments; // K¹t phi (poziomy)
            float phi2 = (j + 1) * 2 * PI / segments;

            // Wierzcho³ki trójk¹ta 1
            float x1 = radius * sin(theta1) * cos(phi1);
            float y1 = radius * sin(theta1) * sin(phi1);
            float z1 = radius * cos(theta1);

            float x2 = radius * sin(theta1) * cos(phi2);
            float y2 = radius * sin(theta1) * sin(phi2);
            float z2 = radius * cos(theta1);

            float x3 = radius * sin(theta2) * cos(phi1);
            float y3 = radius * sin(theta2) * sin(phi1);
            float z3 = radius * cos(theta2);

            // Wierzcho³ki trójk¹ta 2
            float x4 = radius * sin(theta2) * cos(phi2);
            float y4 = radius * sin(theta2) * sin(phi2);
            float z4 = radius * cos(theta2);

            // Rysowanie wype³nionych trójk¹tów (pomarañczowych)
            glColor3f(1.0f, 0.5f, 0.0f); // Pomarañczowy kolor
            glBegin(GL_TRIANGLES);
            glVertex3f(x1, y1, z1);
            glVertex3f(x2, y2, z2);
            glVertex3f(x3, y3, z3);

            glVertex3f(x2, y2, z2);
            glVertex3f(x4, y4, z4);
            glVertex3f(x3, y3, z3);
            glEnd();

            // Rysowanie krawêdzi trójk¹tów (czarnych)
            glColor3f(0.0f, 0.0f, 0.0f); // Czarny kolor
            glBegin(GL_LINE_LOOP);
            glVertex3f(x1, y1, z1);
            glVertex3f(x2, y2, z2);
            glVertex3f(x3, y3, z3);
            glEnd();

            glBegin(GL_LINE_LOOP);
            glVertex3f(x2, y2, z2);
            glVertex3f(x4, y4, z4);
            glVertex3f(x3, y3, z3);
            glEnd();
        }
    }
}
void drawIcosahedron(float radius) {
    const float PI = 3.14159265358979323846f;
    const float H_ANGLE = PI / 180 * 72; // 72 stopni w radianach
    const float V_ANGLE = atanf(1.0f / 2); // K¹t wierzcho³ka

    // Wierzcho³ki dwudziestoœcianu
    float vertices[12][3];

    // Górny wierzcho³ek
    vertices[0][0] = 0;
    vertices[0][1] = radius;
    vertices[0][2] = 0;

    // Dolny wierzcho³ek
    vertices[11][0] = 0;
    vertices[11][1] = -radius;
    vertices[11][2] = 0;

    // 10 wierzcho³ków na œrodku
    for (int i = 1; i <= 5; ++i) {
        vertices[i][0] = radius * cos(V_ANGLE) * cos(i * H_ANGLE);
        vertices[i][1] = radius * sin(V_ANGLE);
        vertices[i][2] = radius * cos(V_ANGLE) * sin(i * H_ANGLE);

        vertices[i + 5][0] = radius * cos(V_ANGLE) * cos((i + 0.5f) * H_ANGLE);
        vertices[i + 5][1] = -radius * sin(V_ANGLE);
        vertices[i + 5][2] = radius * cos(V_ANGLE) * sin((i + 0.5f) * H_ANGLE);
    }

    // Indeksy trójk¹tów
    int indices[20][3] = {
        {0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 5, 1},
        {11, 6, 7}, {11, 7, 8}, {11, 8, 9}, {11, 9, 10}, {11, 10, 6},
        {1, 2, 6}, {2, 3, 7}, {3, 4, 8}, {4, 5, 9}, {5, 1, 10},
        {6, 7, 2}, {7, 8, 3}, {8, 9, 4}, {9, 10, 5}, {10, 6, 1}
    };

    // Rysowanie œcian (pomarañczowe)
    glColor4f(1.0f, 0.5f, 0.0f, 0.5f); // Pomarañczowy kolor
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 3; ++j) {
            glVertex3f(vertices[indices[i][j]][0], vertices[indices[i][j]][1], vertices[indices[i][j]][2]);
        }
    }
    glEnd();

    // Rysowanie krawêdzi (czarne)
    glColor4f(0.0f, 0.0f, 0.0f, 0.5f); // Czarny kolor
    glBegin(GL_LINES);
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 3; ++j) {
            int current = indices[i][j];
            int next = indices[i][(j + 1) % 3]; // Nastêpny wierzcho³ek w trójk¹cie
            glVertex3f(vertices[current][0], vertices[current][1], vertices[current][2]);
            glVertex3f(vertices[next][0], vertices[next][1], vertices[next][2]);
        }
    }
    glEnd();
}