#include <iostream>
#include <GLFW/glfw3.h>
#include <gl/GL.h>
#include <vector>
#include <fstream>
#include <strstream>

#include <algorithm>
#include <list>
#include <Windows.h>

using namespace std;

int width = 1920;
int height = 1080;

struct mat4x4
{
    float m[4][4] = { 0 };
};

struct vec3d 
{
    float x = 0;
    float y = 0;
    float z = 0;
    float w = 1;
};

struct triangle
{
    vec3d p[3];

    float col = 0.0f;
};

int trisLoaded;
int trisDrawn;

struct mesh
{
    vector<triangle> tris;

    bool LoadFromObjFile(string sFilename)
    {
        ifstream f(sFilename);
        if (!f.is_open())
            return false;

        // Local cache of verts
        vector<vec3d> verts;

        while (!f.eof())
        {
            char line[128];
            f.getline(line, 128);

            strstream s;
            s << line;

            char junk;

            if (line[0] == 'v')
            {
                vec3d v;
                s >> junk >> v.x >> v.y >> v.z;
                verts.push_back(v);
            }

            if (line[0] == 'f')
            {
                int f[3];
                s >> junk >> f[0] >> f[1] >> f[2];
                tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });

                trisLoaded++;
            }
        }
        cout << trisLoaded;
        return true;
    }
};


mesh meshCube;
mat4x4 matProj;

vec3d vCamera;
vec3d vLookDir;

float fYaw;

float pX = 0;

int fb = 0;
int ud = 0;

string pKey;

vec3d Matrix_MultiplyVector(mat4x4 &m, vec3d &i)
{
    vec3d v;
    v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
    v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
    v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
    v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
    return v;
}

vec3d Vector_Add(vec3d& v1, vec3d& v2)
{
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

vec3d Vector_Sub(vec3d& v1, vec3d& v2)
{
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

vec3d Vector_Mul(vec3d& v1, float k)
{
    return { v1.x * k, v1.y * k, v1.z * k };
}

vec3d Vector_Div(vec3d& v1, float k)
{
    return { v1.x / k, v1.y / k, v1.z / k };
}

float Vector_DotProduct(vec3d& v1, vec3d& v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

float Vector_Length(vec3d& v)
{
    return sqrtf(Vector_DotProduct(v, v));
}

vec3d Vector_Normalise(vec3d& v)
{
    float l = Vector_Length(v);
    return { v.x / l, v.y / l, v.z / l };
}

vec3d Vector_CrossProduct(vec3d& v1, vec3d& v2)
{
    vec3d v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}

mat4x4 Matrix_MakeIdentity()
{
    mat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 Matrix_MakeRotationZ(float fAngleRad)
{
    mat4x4 matrix;
    matrix.m[0][0] = cosf(fAngleRad);
    matrix.m[0][1] = sinf(fAngleRad);
    matrix.m[1][0] = -sinf(fAngleRad);
    matrix.m[1][1] = cosf(fAngleRad);
    matrix.m[2][2] = 1;
    matrix.m[3][3] = 1;
    return matrix;
}

mat4x4 Matrix_MakeRotationY(float fAngleRad)
{
    mat4x4 matrix;
    matrix.m[0][0] = cosf(fAngleRad);
    matrix.m[0][2] = sinf(fAngleRad);
    matrix.m[1][1] = 1;
    matrix.m[2][0] = -sinf(fAngleRad);
    matrix.m[2][2] = cosf(fAngleRad);
    matrix.m[3][3] = 1;
    return matrix;
}

mat4x4 Matrix_MakeRotationX(float fAngleRad)
{
    mat4x4 matrix;
    matrix.m[0][0] = 1;
    matrix.m[1][1] = cosf(fAngleRad);
    matrix.m[1][2] = sinf(fAngleRad);
    matrix.m[2][1] = -sinf(fAngleRad);
    matrix.m[2][2] = cosf(fAngleRad);
    matrix.m[3][3] = 1;
    return matrix;
}

mat4x4 Matrix_MakeTranslation(float x, float y, float z)
{
    mat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    matrix.m[3][0] = x;
    matrix.m[3][1] = y;
    matrix.m[3][2] = z;
    return matrix;
}

mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
{
    float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
    mat4x4 matrix;
    matrix.m[0][0] = fAspectRatio * fFovRad;
    matrix.m[1][1] = fFovRad;
    matrix.m[2][2] = fFar / (fFar - fNear);
    matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
    matrix.m[2][3] = 1.0f;
    matrix.m[3][3] = 0.0f;
    return matrix;
}

mat4x4 Matrix_MultiplyMatrix(mat4x4& m1, mat4x4& m2)
{
    mat4x4 matrix;
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
    return matrix;
}

mat4x4 Matrix_PointAt(vec3d& pos, vec3d& target, vec3d& up)
{
    vec3d newForward = Vector_Sub(target, pos);
    newForward = Vector_Normalise(newForward);

    vec3d a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));
    vec3d newUp = Vector_Sub(up, a);
    newUp = Vector_Normalise(newUp);

    vec3d newRight = Vector_CrossProduct(newUp, newForward);

    mat4x4 matrix;
    matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 Matrix_QuickInverse(mat4x4& m)
{
    mat4x4 matrix;
    matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
    matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
    matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

vec3d Vector_IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd)
{
    plane_n = Vector_Normalise(plane_n);
    float plane_d = -Vector_DotProduct(plane_n, plane_p);
    float ad = Vector_DotProduct(lineStart, plane_n);
    float bd = Vector_DotProduct(lineEnd, plane_n);
    float t = (-plane_d - ad) / (bd - ad);
    vec3d lineStartToEnd = Vector_Sub(lineEnd, lineStart);
    vec3d lineToIntersect = Vector_Mul(lineStartToEnd, t);
    return Vector_Add(lineStart, lineToIntersect);
}

int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
{
    plane_n = Vector_Normalise(plane_n);

    auto dist = [&](vec3d& p)
    {
        vec3d n = Vector_Normalise(p);
        return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
    };

    vec3d* inside_points[3];  int nInsidePointCount = 0;
    vec3d* outside_points[3]; int nOutsidePointCount = 0;

    float d0 = dist(in_tri.p[0]);
    float d1 = dist(in_tri.p[1]);
    float d2 = dist(in_tri.p[2]);

    if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
    if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
    if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

    if (nInsidePointCount == 0)
    {
        return 0; 
    }

    if (nInsidePointCount == 3)
    {
        out_tri1 = in_tri;

        return 1; 
    }

    if (nInsidePointCount == 1 && nOutsidePointCount == 2)
    {
        out_tri1.col = in_tri.col;

        out_tri1.p[0] = *inside_points[0];

        out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
        out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

        return 1;
    }

    if (nInsidePointCount == 2 && nOutsidePointCount == 1)
    {
        out_tri1.col = in_tri.col;

        out_tri2.col = in_tri.col;

        out_tri1.p[0] = *inside_points[0];
        out_tri1.p[1] = *inside_points[1];
        out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

        out_tri2.p[0] = *inside_points[1];
        out_tri2.p[1] = out_tri1.p[2];
        out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

        return 2;
    }
}

float DegToRad(float a)
{
    return { 180.0f / 3.14159f * a };
}

void createTris()
{
    meshCube.LoadFromObjFile("terraintest.obj");

    matProj = Matrix_MakeProjection(90.0f, (float)height / (float)width, 0.1f, 1000.0f);
}

void drawTriangle(float x1, float y1, float x2, float y2, float x3, float y3, float dp)
{
    trisDrawn++;
    glColor3f(0.5, 0.5, dp);
    glBegin(GL_TRIANGLES); glVertex2f(x1, y1); glVertex2f(x2, y2); glVertex2f(x3, y3); glEnd();
    //glColor3f(1.0f, 1.0f, 1.0f);
    //glLineWidth(2); glBegin(GL_LINES); glVertex2f(x1, y1); glVertex2f(x2, y2); glEnd();
    //glLineWidth(2); glBegin(GL_LINES); glVertex2f(x2, y2); glVertex2f(x3, y3); glEnd();
    //glLineWidth(2); glBegin(GL_LINES); glVertex2f(x3, y3); glVertex2f(x1, y1); glEnd();
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_UP)
        switch (action)
            {
            case GLFW_PRESS:
                ud = 2;
                break;

            case GLFW_REPEAT:
                ud = 2;
                break;

            case GLFW_RELEASE:
                ud = 0;
                break;
        }
    if (key == GLFW_KEY_DOWN)
        switch (action)
            {
            case GLFW_PRESS:
                ud = 1;
                break;

            case GLFW_REPEAT:
                ud = 1;
                break;

            case GLFW_RELEASE:
                ud = 0;
                break;
        }
    if (key == GLFW_KEY_D)
        switch (action)
        {
            case GLFW_PRESS:
                fYaw += 0.01;
                break;

            case GLFW_REPEAT:
                fYaw += 0.01;
                break;
        }

    if (key == GLFW_KEY_A)
        switch (action)
        {
        case GLFW_PRESS:
            fYaw -= 0.01;
            break;

        case GLFW_REPEAT:
            fYaw -= 0.01;
            break;
        }
    if (key == GLFW_KEY_W)
    {
        switch (action)
        {
            case GLFW_PRESS:
                fb = 2;
                break;

            case GLFW_REPEAT:
                fb = 2;
                break;
            
            case GLFW_RELEASE:
                fb = 0;
                break;
        }
    }
    else if (key == GLFW_KEY_S)
        switch (action)
        {
        case GLFW_PRESS:
            fb = 1;
            break;

        case GLFW_REPEAT:
            fb = 1;
            break;

        case GLFW_RELEASE:
            fb = 0;
            break;
        }
}

int main(void)
{
    GLFWwindow* window;

    if (!glfwInit()) {
        return -1;
    }

    window = glfwCreateWindow(width, height, "Testi", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);

    glClearColor(0.1f, 0.1f, 0.1f, 0);
    glOrtho(0, width, height, 0, 0, 1);

    createTris();

    float zAngle = DegToRad(0.0f);
    float xAngle = DegToRad(0.0f);

    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        //glLoadIdentity();

        vec3d vForward = Vector_Mul(vLookDir, 0.1f);

        if (fb == 2)
            vCamera = Vector_Add(vCamera, vForward);
        else if (fb == 1)
            vCamera = Vector_Sub(vCamera, vForward);

        if (ud == 2)
            vCamera.y += 1;
        else if (ud == 1)
            vCamera.y -= 1;

        mat4x4 matRotZ, matRotX;

        //vCamera.x = pX;

        //zAngle += 0.005f;
        //xAngle += 0.001f;

        matRotZ = Matrix_MakeRotationZ(zAngle);
        matRotX = Matrix_MakeRotationX(xAngle);

        mat4x4 matTrans;
        matTrans = Matrix_MakeTranslation(0.0f, -8.0f, 8.0f);

        mat4x4 matWorld;
        matWorld = Matrix_MakeIdentity();
        matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);
        matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);

        vec3d vUp = { 0,1,0 };
        vec3d vTarget = { 0,0,1 };
        mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
        vLookDir = Matrix_MultiplyVector(matCameraRot, vTarget);
        vTarget = Vector_Add(vCamera, vLookDir);
        mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);

        mat4x4 matView = Matrix_QuickInverse(matCamera);

        vector<triangle> vecTrianglesToRaster;

        for (auto tri : meshCube.tris)
        {
            triangle triProjected, triTransformed, triViewed;

            triTransformed.p[0] = Matrix_MultiplyVector(matWorld, tri.p[0]);
            triTransformed.p[1] = Matrix_MultiplyVector(matWorld, tri.p[1]);
            triTransformed.p[2] = Matrix_MultiplyVector(matWorld, tri.p[2]);

            vec3d normal, line1, line2;

            line1 = Vector_Sub(triTransformed.p[1], triTransformed.p[0]);
            line2 = Vector_Sub(triTransformed.p[2], triTransformed.p[0]);

            normal = Vector_CrossProduct(line1, line2);

            normal = Vector_Normalise(normal);

            float l = sqrtf(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
            normal.x /= l; normal.y /= l; normal.z /= l;

            vec3d vCameraRay = Vector_Sub(triTransformed.p[0], vCamera);

            if (Vector_DotProduct(normal, vCameraRay) < 0.0f)
            {
                vec3d light_direction = { 0.0f, 1.0f, -1.0f };
                light_direction = Vector_Normalise(light_direction);

                float dp = normal.x * light_direction.x + normal.y * light_direction.y + normal.z * light_direction.z;

                triViewed.p[0] = Matrix_MultiplyVector(matView, triTransformed.p[0]);
                triViewed.p[1] = Matrix_MultiplyVector(matView, triTransformed.p[1]);
                triViewed.p[2] = Matrix_MultiplyVector(matView, triTransformed.p[2]);

                int nClippedTriangles = 0;
                triangle clipped[2];
                nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);

                for (int n = 0; n < nClippedTriangles; n++)
                {


                    triProjected.p[0] = Matrix_MultiplyVector(matProj, clipped[n].p[0]);
                    triProjected.p[1] = Matrix_MultiplyVector(matProj, clipped[n].p[1]);
                    triProjected.p[2] = Matrix_MultiplyVector(matProj, clipped[n].p[2]);
                    triProjected.col = dp;

                    triProjected.p[0] = Vector_Div(triProjected.p[0], triProjected.p[0].w);
                    triProjected.p[1] = Vector_Div(triProjected.p[1], triProjected.p[1].w);
                    triProjected.p[2] = Vector_Div(triProjected.p[2], triProjected.p[2].w);

                    triProjected.p[0].x *= -1.0f;
                    triProjected.p[1].x *= -1.0f;
                    triProjected.p[2].x *= -1.0f;
                    triProjected.p[0].y *= -1.0f;
                    triProjected.p[1].y *= -1.0f;
                    triProjected.p[2].y *= -1.0f;

                    // Scaling
                    vec3d vOffsetView = { 1,1,0 };
                    triProjected.p[0] = Vector_Add(triProjected.p[0], vOffsetView);
                    triProjected.p[1] = Vector_Add(triProjected.p[1], vOffsetView);
                    triProjected.p[2] = Vector_Add(triProjected.p[2], vOffsetView);

                    triProjected.p[0].x *= 0.5f * (float)width;
                    triProjected.p[0].y *= 0.5f * (float)height;
                    triProjected.p[1].x *= 0.5f * (float)width;
                    triProjected.p[1].y *= 0.5f * (float)height;
                    triProjected.p[2].x *= 0.5f * (float)width;
                    triProjected.p[2].y *= 0.5f * (float)height;

                    vecTrianglesToRaster.push_back(triProjected);

                }
            }
        }
        sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(), [](triangle& t1, triangle& t2)
        {
            float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
            float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
            return z1 > z2;
        });

        for (auto& triToRaster : vecTrianglesToRaster)
            {
                triangle clipped[2];
                list<triangle> listTriangles;

                listTriangles.push_back(triToRaster);
                int nNewTriangles = 1;

                for (int p = 0; p < 4; p++)
                {
                    int nTrisToAdd = 0;
                    while (nNewTriangles > 0)
                    {
                        triangle test = listTriangles.front();
                        listTriangles.pop_front();
                        nNewTriangles--;

                        switch (p)
                        {
                        case 0:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                        case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)height - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                        case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                        case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)width - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                        }

                        for (int w = 0; w < nTrisToAdd; w++)
                            listTriangles.push_back(clipped[w]);
                    }
                    nNewTriangles = listTriangles.size();
                }
                for (auto& t : listTriangles)
                {
                    drawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, t.col);
                }
            }


        cout << "Tris drawn: " << trisDrawn;

        

        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

