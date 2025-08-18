
#include <GLFW/glfw3.h>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cmath>
#define M_PI 3.1415

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void renderScene();
void drawCuboid(float width, float height, float depth);
void drawOrangeSphereWithBlackEdges(float radius, int segments);

struct node {
    glm::vec3 position = glm::vec3(0, 0, 0); // Domyslnie ustawione na (0,0,0)
    glm::vec3 velocity = glm::vec3(0, 0, 0);
};
std::vector<node> nodes;

struct Cuboid_dimensions {
    float width;
    float height;
    float depth;
};

struct Triangle {
    int indices[3];
};
std::vector<Triangle> triangles;

void drawIcosahedron(float radius, std::vector<node>, std::vector<Triangle>);

void updatePhysics(float dt, float cubeWidth, float cubeHeight, float cubeDepth) {
    const float cuboidHalfWidth = cubeWidth / 2.0f;
    const float cuboidHalfHeight = cubeHeight / 2.0f;
    const float cuboidHalfDepth = cubeDepth / 2.0f;
    const float elasticity = 1.0f; // wspoolczynnik sprezystosci (0.8 = 80% energii zachowane)

    for (auto& node : nodes) {
        // Aktualizacja pozycji
        node.position = node.position + node.velocity * dt;

        // Sprawdzanie kolizji ze scianami prostopadloscianu i odbicia
        if (node.position.x <= -cuboidHalfWidth || node.position.x >= cuboidHalfWidth) {
            node.velocity.x = -node.velocity.x * elasticity;
            // Korekta pozycji aby nie utknac w scianie
            node.position.x = (node.position.x < 0) ? -cuboidHalfWidth + 0.001f : cuboidHalfWidth - 0.001f;
        }

        if (node.position.y <= -cuboidHalfHeight || node.position.y >= cuboidHalfHeight) {
            node.velocity.y = -node.velocity.y * elasticity;
            node.position.y = (node.position.y < 0) ? -cuboidHalfHeight + 0.001f : cuboidHalfHeight - 0.001f;
        }

        if (node.position.z <= -cuboidHalfDepth || node.position.z >= cuboidHalfDepth) {
            node.velocity.z = -node.velocity.z * elasticity;
            node.position.z = (node.position.z < 0) ? -cuboidHalfDepth + 0.001f : cuboidHalfDepth - 0.001f;
        }
    }
}

float calculateTriangleArea(const std::vector<node>& nodes, int a, int b, int c) {
    glm::vec3 ab = nodes[b].position - nodes[a].position;
    glm::vec3 ac = nodes[c].position - nodes[a].position;
    glm::vec3 cross = glm::cross(ab, ac);
    return 0.5f * glm::length(cross);
}

int addMidpoint(std::vector<node>& nodes, int a, int b) {
    node midpoint;
    midpoint.position = (nodes[a].position + nodes[b].position) * 0.5f;
    midpoint.velocity = (nodes[a].velocity + nodes[b].velocity) * 0.5f;
    nodes.push_back(midpoint);
    return (int)(nodes.size() - 1);
}

void refineIcosahedron(std::vector<node>& nodes, std::vector<Triangle>& triangles, float maxArea) {
    std::vector<Triangle> newTriangles;

    for (const auto& tri : triangles) {
        int a = tri.indices[0], b = tri.indices[1], c = tri.indices[2];
        float area = calculateTriangleArea(nodes, a, b, c);

        if (area > maxArea) {
            // Dodaj punkty srodkowe krawedzi
            int ab = addMidpoint(nodes, a, b);
            int bc = addMidpoint(nodes, b, c);
            int ca = addMidpoint(nodes, c, a);

            // Podziel trojkat na 4 mniejsze
            newTriangles.push_back({ {a, ab, ca} });
            newTriangles.push_back({ {b, bc, ab} });
            newTriangles.push_back({ {c, ca, bc} });
            newTriangles.push_back({ {ab, bc, ca} });
        }
        else {
            newTriangles.push_back(tri);
        }
    }

    triangles = newTriangles;
}

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
bool first = true;

float radius = 0.5f;
const float H_ANGLE = M_PI / 180 * 72; // 72 stopni w radianach
const float V_ANGLE = atanf(1.0f / 2); // Kat wierzcholka

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
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

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
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    float dt = 0.01;
    static std::vector<Triangle> triangles;
    if (first) {
        node buf;
        buf.position = glm::vec3(0.0f, radius, 0.0f);
        nodes.push_back(buf);

        for (int i = 1; i <= 5; ++i) {
            buf.position = glm::vec3(radius * cos(V_ANGLE) * cos(i * H_ANGLE), radius * sin(V_ANGLE), radius * cos(V_ANGLE) * sin(i * H_ANGLE));
            nodes.push_back(buf);
        }
        for (int i = 6; i <= 10; ++i) {
            buf.position = glm::vec3(radius * cos(V_ANGLE) * cos((i + 0.5f) * H_ANGLE), -radius * sin(V_ANGLE), radius * cos(V_ANGLE) * sin((i + 0.5f) * H_ANGLE));
            nodes.push_back(buf);
        }

        // Dolny wierzcholek
        buf.position = glm::vec3(0.0f, -radius, 0.0f);
        nodes.push_back(buf);
        first = false;

        for (int i = 0; i < nodes.size(); ++i) {
            nodes[i].velocity = nodes[i].position * 0.1f;
        }

        const int initTris[20][3] = {
            {0,1,2}, {0,2,3}, {0,3,4}, {0,4,5}, {0,5,1},
            {11,6,7}, {11,7,8}, {11,8,9}, {11,9,10}, {11,10,6},
            {1,2,6}, {2,3,7}, {3,4,8}, {4,5,9}, {5,1,10},
            {6,7,2}, {7,8,3}, {8,9,4}, {9,10,5}, {10,6,1}
        };

        for (int i = 0; i < 20; ++i) {
            triangles.push_back({ {initTris[i][0], initTris[i][1], initTris[i][2]} });
        }
    }

    static int frameCount = 0;
    if (frameCount % 10 == 0) {
        refineIcosahedron(nodes, triangles, 0.2f);
    }
    frameCount++;

    //updatePhysics(dt);
    Cuboid_dimensions Cube;
    Cube.width = 8.0f;
    Cube.height = 6.0f;
    Cube.depth = 10.0f;
    updatePhysics(dt, Cube.width, Cube.height, Cube.depth);

    drawIcosahedron(0.5f, nodes, triangles);

    // Basen
    glColor4f(0.0f, 0.0f, 1.0f, 0.1f);
    drawCuboid(Cube.width, Cube.height, Cube.depth);
}

void drawCuboid(float width, float height, float depth) {
    float halfWidth = width / 2.0f;
    float halfHeight = height / 2.0f;
    float halfDepth = depth / 2.0f;

    // Definiowanie wierzcholkow prostopadloscianu
    glm::vec3 vertices[] = {
        // Przednia sciana
        { -halfWidth, -halfHeight,  halfDepth },
        {  halfWidth, -halfHeight,  halfDepth },
        {  halfWidth,  halfHeight,  halfDepth },
        { -halfWidth,  halfHeight,  halfDepth },

        // Tylna sciana
        { -halfWidth, -halfHeight, -halfDepth },
        {  halfWidth, -halfHeight, -halfDepth },
        {  halfWidth,  halfHeight, -halfDepth },
        { -halfWidth,  halfHeight, -halfDepth }
    };

    // Indeksy wierzcholkow dla scian (czworokatow)
    int faces[6][4] = {
        {0, 1, 2, 3}, // Przednia sciana
        {4, 5, 6, 7}, // Tylna sciana
        {0, 3, 7, 4}, // Lewa sciana
        {1, 2, 6, 5}, // Prawa sciana
        {0, 1, 5, 4}, // Dolna sciana
        {2, 3, 7, 6}  // Gorna sciana
    };

    // Rysowanie scian (czworokatow)
    glBegin(GL_QUADS);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            glVertex3f(vertices[faces[i][j]].x, vertices[faces[i][j]].y, vertices[faces[i][j]].z);
        }
    }
    glEnd();

    // Rysowanie krawedzi (linii)
    glColor3f(0.0f, 0.0f, 0.0f);
    glBegin(GL_LINES);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            int current = faces[i][j];
            int next = faces[i][(j + 1) % 4];
            glVertex3f(vertices[current].x, vertices[current].y, vertices[current].z);
            glVertex3f(vertices[next].x, vertices[next].y, vertices[next].z);
        }
    }
    glEnd();
}

void drawOrangeSphereWithBlackEdges(float radius, int segments) {
    const float PI = 3.14159265358979323846f;

    for (int i = 0; i < segments; ++i) {
        float theta1 = i * PI / segments;
        float theta2 = (i + 1) * PI / segments;

        for (int j = 0; j < segments; ++j) {
            float phi1 = j * 2 * PI / segments;
            float phi2 = (j + 1) * 2 * PI / segments;

            // Wierzcholki trojkata 1
            float x1 = radius * sin(theta1) * cos(phi1);
            float y1 = radius * sin(theta1) * sin(phi1);
            float z1 = radius * cos(theta1);

            float x2 = radius * sin(theta1) * cos(phi2);
            float y2 = radius * sin(theta1) * sin(phi2);
            float z2 = radius * cos(theta1);

            float x3 = radius * sin(theta2) * cos(phi1);
            float y3 = radius * sin(theta2) * sin(phi1);
            float z3 = radius * cos(theta2);

            // Wierzcholki trojkata 2
            float x4 = radius * sin(theta2) * cos(phi2);
            float y4 = radius * sin(theta2) * sin(phi2);
            float z4 = radius * cos(theta2);

            // Rysowanie wypelnionych trojkatow
            glColor3f(1.0f, 0.5f, 0.0f);
            glBegin(GL_TRIANGLES);
            glVertex3f(x1, y1, z1);
            glVertex3f(x2, y2, z2);
            glVertex3f(x3, y3, z3);

            glVertex3f(x2, y2, z2);
            glVertex3f(x4, y4, z4);
            glVertex3f(x3, y3, z3);
            glEnd();

            // Rysowanie krawedzi trojkatow
            glColor3f(0.0f, 0.0f, 0.0f);
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

void drawIcosahedron(float radius, std::vector<node> nodes, std::vector<Triangle> triangles) {
    // Rysowanie scian (pomarañczowe)
    glColor4f(1.0f, 0.5f, 0.0f, 0.5f);
    glBegin(GL_TRIANGLES);
    for (const auto& tri : triangles) {
        for (int j = 0; j < 3; ++j) {
            glVertex3f(nodes[tri.indices[j]].position.x,
                nodes[tri.indices[j]].position.y,
                nodes[tri.indices[j]].position.z);
        }
    }
    glEnd();

    // Rysowanie krawedzi (czarne)
    glColor4f(0.0f, 0.0f, 0.0f, 0.5f);
    glBegin(GL_LINES);
    for (const auto& tri : triangles) {
        for (int j = 0; j < 3; ++j) {
            int current = tri.indices[j];
            int next = tri.indices[(j + 1) % 3];
            glVertex3f(nodes[current].position.x, nodes[current].position.y, nodes[current].position.z);
            glVertex3f(nodes[next].position.x, nodes[next].position.y, nodes[next].position.z);
        }
    }
    glEnd();
}