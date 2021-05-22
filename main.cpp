#include <iostream>
#include <stack>
#include <vector>
#include <sstream>
#include <cmath>
#include <ctime>

#include "bitmap_image.hpp"

using namespace std;

#define N 4
#define PI (2*acos(0.0))

class Vector{
public:
    double x, y, z;

    Vector() : x(0), y(0), z(0) {}

    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    explicit Vector(const string& line) : x(0), y(0), z(0) {
        stringstream ss(line);
        ss >> x;
        ss >> y;
        ss >> z;
    }

    void normalize() {
        double value = sqrt(x*x + y*y + z*z);
        x = x/value;
        y = y/value;
        z = z/value;
    }

    double dotProduct(Vector* v) const {
        return x*v->x + y*v->y + z*v->z;
    }

    Vector* crossProduct(Vector* v) const {
        double a = y * v->z - z * v->y;
        double b = z * v->x - x * v->z;
        double c = x * v->y - y * v->x;

        return new Vector(a, b, c);
    }

    Vector* operator + (const Vector &v) const {
        double a = x + v.x;
        double b = y + v.y;
        double c = z + v.z;

        return new Vector(a, b, c);
    }

    Vector* operator - (const Vector &v) const {
        double a = x - v.x;
        double b = y - v.y;
        double c = z - v.z;

        return new Vector(a, b, c);
    }

    string toString() const {
        return to_string(x) + " " + to_string(y) + " " + to_string(z);
    }
};

class Color{
public:
    int r, g, b;

    Color()  : r(0), g(0), b(0) {}

    Color(int r, int g, int b) : r(r), g(g), b(b) {}
};

class Triangle{
public:
    Vector *a, *b, *c;
    Color color;

    Triangle(Vector *a, Vector *b, Vector *c) : a(a), b(b), c(c) {}

    Triangle(Vector *a, Vector *b, Vector *c, const Color &color) : a(a), b(b), c(c), color(color) {}

    string toString() const {
        return a->toString() + "\n" + b->toString() + "\n" + c->toString();
    }
};

double** initializeMatrix() {
    double** matrix;
    matrix = new double*[N];

    for (int i = 0; i < N; ++i) {
        matrix[i] = new double[N];
        for (int j = 0; j < N; ++j) {
            matrix[i][j] = (i==j) ? 1 : 0;
        }
    }

    return matrix;
}

double** product(double** a, double** b) {
    double** result;
    result = new double*[N];

    for (int i = 0; i < N; ++i) {
        result[i] = new double[N];
        for (int j = 0; j < N; ++j) {
            result[i][j] = 0;
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                result[i][j] += a[i][k]*b[k][j];
            }
        }
    }

    return result;
}

Vector* transform(double** matrix, const Vector* p) {
    double resultMatrix[] = {0, 0, 0, 0} ;
    double pointMatrix[] = {p->x, p->y, p->z, 1};

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            resultMatrix[i] += matrix[i][j]*pointMatrix[j];
        }
    }

    return new Vector(resultMatrix[0]/resultMatrix[3], resultMatrix[1]/resultMatrix[3], resultMatrix[2]/resultMatrix[3]);
}

Vector* rotate(Vector* v, Vector* axis, double angle) {
    double theta = angle * (PI/180);
    double dotProduct = axis->dotProduct(v);
    Vector* crossProduct = axis->crossProduct(v);

    double a = v->x * cos(theta) + crossProduct->x * sin(theta) + axis->x * dotProduct * (1 - cos(theta));
    double b = v->y * cos(theta) + crossProduct->y * sin(theta) + axis->y * dotProduct * (1 - cos(theta));
    double c = v->z * cos(theta) + crossProduct->z * sin(theta) + axis->z * dotProduct * (1 - cos(theta));

    return new Vector(a, b, c);
}

double** makeTranslateMatrix(double tx, double ty, double tz) {
    double** translateMatrix = initializeMatrix();
    translateMatrix[0][3] = tx;
    translateMatrix[1][3] = ty;
    translateMatrix[2][3] = tz;

    return translateMatrix;
}

double** makeScaleMatrix(double sx, double sy, double sz) {
    double** scaleMatrix = initializeMatrix();
    scaleMatrix[0][0] = sx;
    scaleMatrix[1][1] = sy;
    scaleMatrix[2][2] = sz;

    return scaleMatrix;
}

double** makeRotationMatrix(Vector* axis, double angle) {
    Vector* c1 = rotate(new Vector(1, 0, 0), axis, angle);
    Vector* c2 = rotate(new Vector(0, 1, 0), axis, angle);
    Vector* c3 = rotate(new Vector(0, 0, 1), axis, angle);

    double** result = initializeMatrix();
    result[0][0] = c1->x;
    result[1][0] = c1->y;
    result[2][0] = c1->z;

    result[0][1] = c2->x;
    result[1][1] = c2->y;
    result[2][1] = c2->z;

    result[0][2] = c3->x;
    result[1][2] = c3->y;
    result[2][2] = c3->z;

    return result;
}

void printMatrix(double** matrix, const char *message = nullptr) {
    if(message != nullptr) cout << message << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

Triangle* transformTriangle(Triangle* triangle, double** transformMatrix) {
    return new Triangle(transform(transformMatrix, triangle->a), transform(transformMatrix, triangle->b), transform(transformMatrix, triangle->c));
}

double** getViewMatrix(Vector* eye, Vector* look, Vector* up) {
    Vector *l, *r, *u;
    l = *look-*eye;
    l->normalize();

    r = l->crossProduct(up);
    r->normalize();

    u = r->crossProduct(l);

    double** eyeTranslateMatrix = initializeMatrix();
    eyeTranslateMatrix[0][3] = -eye->x;
    eyeTranslateMatrix[1][3] = -eye->y;
    eyeTranslateMatrix[2][3] = -eye->z;

    double** eyeRotateMatrix = initializeMatrix();
    eyeRotateMatrix[0][0] = r->x;
    eyeRotateMatrix[0][1] = r->y;
    eyeRotateMatrix[0][2] = r->z;

    eyeRotateMatrix[1][0] = u->x;
    eyeRotateMatrix[1][1] = u->y;
    eyeRotateMatrix[1][2] = u->z;

    eyeRotateMatrix[2][0] = -l->x;
    eyeRotateMatrix[2][1] = -l->y;
    eyeRotateMatrix[2][2] = -l->z;

    return product(eyeRotateMatrix, eyeTranslateMatrix);
}

double** getProjectionMatrix(double fovY, double aspectRatio, double near, double far, double** viewMatrix) {
    double fovX = fovY * aspectRatio;
    double _t = near * tan((fovY/2)*(PI/180));
    double _r = near * tan((fovX/2)*(PI/180));

    double** projectionMatrix = initializeMatrix();

    projectionMatrix[0][0] = near/_r;
    projectionMatrix[1][1] = near/_t;
    projectionMatrix[2][2] = -(far+near)/(far-near);
    projectionMatrix[2][3] = -(2*far*near)/(far-near);
    projectionMatrix[3][2] = -1;
    projectionMatrix[3][3] = 0;

    return product(projectionMatrix, viewMatrix);
}

int main() {
    srand(time(nullptr));
    Vector *eye, *look, *up;
    double fovY, aspectRatio, near, far;

    string line;
    ifstream inputFile;
    ofstream outputFile;

    string dir = "test/1/";
    string sceneFileName = dir+"scene.txt";
    string configFileName = dir+"config.txt";
    string modelOutputFileName = "out_stage1.txt";
    string viewOutputFileName = "out_stage2.txt";
    string projectionOutputFileName = "out_stage3.txt";
    string zBufferOutputFileName = "out_z_buffer.txt";
    string imageFileName = "out.bmp";

    inputFile.open(sceneFileName);

    getline(inputFile, line);
    eye = new Vector(line);

    getline(inputFile, line);
    look = new Vector(line);

    getline(inputFile, line);
    up = new Vector(line);

    inputFile >> fovY;
    inputFile >> aspectRatio;
    inputFile >> near;
    inputFile >> far;

    stack<double**> s;
    stack<int> pushPoints;
    vector<Triangle*> triangles;

    s.push(initializeMatrix());

    while(getline(inputFile, line)) {
        if (line == "triangle") {
            Vector *a, *b, *c;
            getline(inputFile, line);
            a = new Vector(line);

            getline(inputFile, line);
            b = new Vector(line);

            getline(inputFile, line);
            c = new Vector(line);

            triangles.push_back(transformTriangle(new Triangle(a, b, c), s.top()));
        }
        else if (line == "translate") {
            double tx, ty, tz;
            inputFile >> tx;
            inputFile >> ty;
            inputFile >> tz;

            s.push(product(s.top(), makeTranslateMatrix(tx, ty, tz)));
        }
        else if (line == "rotate") {
            double angle, x, y, z;

            inputFile >> angle;
            inputFile >> x;
            inputFile >> y;
            inputFile >> z;

            auto* axis = new Vector(x, y, z);
            axis->normalize();

            s.push(product(s.top(), makeRotationMatrix(axis, angle)));
        }
        else if (line == "scale") {
            double sx, sy, sz;
            inputFile >> sx;
            inputFile >> sy;
            inputFile >> sz;

            s.push(product(s.top(), makeScaleMatrix(sx, sy, sz)));
        }
        else if (line == "push") {
            pushPoints.push((int)s.size());
            //printMatrix(s.top());
        }
        else if (line == "pop") {
            if(pushPoints.empty()) {
                cout << "Error: pop before push(pop ignored)" << endl;
            } else {
                //printMatrix(s.top());
                int pushPoint = pushPoints.top();
                while(s.size() != pushPoint) {
                    s.pop();
                }
                pushPoints.pop();
                //printMatrix(s.top());
            }
        }
        else if (line == "end") {
            break;
        }

    }
    inputFile.close();

    outputFile.open(modelOutputFileName);
    for (auto & triangle : triangles) {
        outputFile << triangle->toString() << endl << endl;
    }
    outputFile.close();


    //task 2
    double** viewMatrix = getViewMatrix(eye, look, up);
    //printMatrix(viewMatrix, "viewMatrix:");

    outputFile.open(viewOutputFileName);
    for (auto & triangle : triangles) {
        outputFile << transformTriangle(triangle, viewMatrix)->toString() << endl << endl;
    }
    outputFile.close();

    //task 3
    double** projectionMatrix = getProjectionMatrix(fovY, aspectRatio, near, far, viewMatrix);
    //printMatrix(projectionMatrix);

    outputFile.open(projectionOutputFileName);
    for (auto & triangle : triangles) {
        outputFile << transformTriangle(triangle, projectionMatrix)->toString() << endl << endl;
    }
    outputFile.close();


    //task 4
    int screenWidth, screenHeight;
    double x_min, x_max, y_min, y_max, z_min, z_max;
    double dx, dy, top_y, bottom_y, left_x, right_x;

    vector<Triangle*> projectedTriangles;
    projectedTriangles.reserve(triangles.size());

    inputFile.open(projectionOutputFileName);
    while(getline(inputFile, line)) {
        auto* a = new Vector(line);

        getline(inputFile, line);
        auto* b = new Vector(line);

        getline(inputFile, line);
        auto* c = new Vector(line);

        Color color(rand()%255, rand()%255, rand()%255);
        projectedTriangles.push_back(new Triangle(a, b, c, color));
    }
    inputFile.close();

    inputFile.open(configFileName);
    inputFile >> screenWidth;
    inputFile >> screenHeight;

    inputFile >> x_min;
    inputFile >> y_min;
    inputFile >> z_min >> z_max;

    inputFile.close();

    x_max = -x_min;
    y_max = -y_min;

    dx = (x_max - x_min)/screenWidth;
    dy = (y_max - y_min)/screenHeight;

    top_y = y_max - (dy/2);
    bottom_y = y_min + (dy/2);
    left_x = x_min + (dx/2);
    right_x = x_max - (dx/2);

    //cout << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << endl;

    double** z_buffer;
    Color** frame;

    z_buffer = new double*[screenWidth];
    frame = new Color*[screenWidth];
    for (int i = 0; i < screenWidth; ++i) {
        z_buffer[i] = new double[screenHeight];
        frame[i] = new Color[screenHeight];
        for (int j = 0; j < screenHeight; ++j) {
            z_buffer[i][j] = z_max;
            frame[i][j] = Color(255, 255, 255);
        }
    }

    for (auto & triangle : projectedTriangles) {
        //cout << triangle->toString() << endl << endl;
        double max_point_x, max_point_y, min_point_y, min_point_x;
        max_point_x = max(max(triangle->a->x, triangle->b->x), triangle->c->x);
        max_point_y = max(max(triangle->a->y, triangle->b->y), triangle->c->y);
        min_point_x = min(min(triangle->a->x, triangle->b->x), triangle->c->x);
        min_point_y = min(min(triangle->a->y, triangle->b->y), triangle->c->y);

        max_point_x = max_point_x > right_x ? right_x : max_point_x;
        max_point_y = max_point_y > top_y ? top_y : max_point_y;
        min_point_x = min_point_x < left_x ? left_x : min_point_x;
        min_point_y = min_point_y < bottom_y ? bottom_y : min_point_y;

        int top_scanLine, bottom_scanLine;
        int left_scanLine, right_scanLine;

        top_scanLine = round((top_y - max_point_y) / dy);
        bottom_scanLine = round(screenHeight - (min_point_y - bottom_y) / dy);
        left_scanLine = round((min_point_x - left_x) / dy);
        right_scanLine = round(screenWidth - (right_x - max_point_x) / dx);

        /*cout << "max point on Y: " << max_point_y << endl;
        cout << "min point on Y: " << min_point_y << endl;
        cout << "max point on X: " << max_point_x << endl;
        cout << "min point on X: " << min_point_x << endl << endl;

        cout << "top_scanLine: " << top_scanLine << endl;
        cout << "bottom_scanLine: " << bottom_scanLine << endl;
        cout << "left_scanLine: " << left_scanLine << endl;
        cout << "right_scanLine: " << right_scanLine << endl;
        cout << endl;*/

        for (int i = top_scanLine; i <= bottom_scanLine; ++i) {
            double y_s = max_point_y - i*dy;
            double z_a = triangle->a->z - (triangle->a->z - triangle->b->z)*((triangle->a->y - y_s)/(triangle->a->y - triangle->b->y));
            double z_b = triangle->a->z - (triangle->a->z - triangle->c->z)*((triangle->a->y - y_s)/(triangle->a->y - triangle->c->y));
            for (int j = left_scanLine; j <= right_scanLine; ++j) {
                double x_p = min_point_x + j*dx;
                double z_p = z_b - (z_b - z_a)*((max_point_x - x_p)/(max_point_x-min_point_x));

                if(z_p < z_buffer[i][j] && z_p > z_min) {
                    //cout << "i: " << i << " j: " << j << " ";
                    //cout << "z_p: " << z_p << " z_buffer[i][j]: " << z_buffer[i][j] << endl;
                    //cout << "r: " << triangle->color.r << " g: " << triangle->color.g << " b: " << triangle->color.b << endl;
                    z_buffer[i][j] = z_p;
                    frame[i][j] = triangle->color;
                }
            }
        }

    }

    //image save
    bitmap_image image(screenWidth, screenHeight);
    for(int i=0;i<screenWidth;i++){
        for(int j=0;j<screenHeight;j++){
            image.set_pixel(i,j,frame[i][j].r,frame[i][j].g,frame[i][j].b);
        }
    }

    image.save_image(imageFileName);

    //z_buffer save
    outputFile.open(zBufferOutputFileName);
    for(int i=0;i<screenWidth;i++){
        for(int j=0;j<screenHeight;j++){
            if(z_buffer[i][j] < z_max) outputFile << z_buffer[i][j];
            outputFile << " ";
        }
        outputFile << endl ;
    }
    outputFile.close();

    //free memory
    for (int i = 0; i < screenWidth; ++i) {
        delete[] z_buffer[i];
        delete[] frame[i];
    }
    delete[] z_buffer;
    delete[] frame;


    return 0;
}
