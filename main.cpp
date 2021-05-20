#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <sstream>
#include <cmath>

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

    Vector* operator + (Vector* v) const {
        double a = x + v->x;
        double b = y + v->y;
        double c = z + v->z;

        return new Vector(a, b, c);
    }

    Vector* operator - (Vector* v) const {
        double a = x - v->x;
        double b = y - v->y;
        double c = z - v->z;

        return new Vector(a, b, c);
    }

    string toString() const {
        return to_string(x) + " " + to_string(y) + " " + to_string(z);
    }
};

class Triangle{
public:
    Vector *a, *b, *c;
    int* color;

    Triangle(Vector *a, Vector *b, Vector *c) : a(a), b(b), c(c) {
        color = new int[3];
        for (int i = 0; i < 3; ++i) color[i] = 0;
    }

    void setColor(int _r, int _g, int _b) const {
        color[0] = _r;
        color[1] = _g;
        color[2] = _b;
    }

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
    l = *look-eye;
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


class GL{
public:
    double** viewMatrix;
    double** projectionMatrix;
    stack<double**> transformMatrices;
    stack<int> pushPoints;

    void initialize() {
        viewMatrix = initializeMatrix();
        projectionMatrix = initializeMatrix();

        transformMatrices.push(initializeMatrix());
    }

    void setView(Vector* eye, Vector* look, Vector* up) {
        Vector *l, *r, *u;
        l = *look-eye;
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

        viewMatrix = product(eyeRotateMatrix, eyeTranslateMatrix);
    }

    void setProjection(double fovY, double aspectRatio, double near, double far) {
        double fovX = fovY * aspectRatio;
        double _t = near * tan((fovY/2)*(PI/180));
        double _r = near * tan((fovX/2)*(PI/180));

        projectionMatrix = initializeMatrix();

        projectionMatrix[0][0] = near/_r;
        projectionMatrix[1][1] = near/_t;
        projectionMatrix[2][2] = -(far+near)/(far-near);
        projectionMatrix[2][3] = -(2*far*near)/(far-near);
        projectionMatrix[3][2] = -1;
        projectionMatrix[3][3] = 0;

        projectionMatrix = product(projectionMatrix, viewMatrix);
    }


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
    }

    Vector* transform(const Vector* p, double** matrix) {
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

    void push() {
        pushPoints.push((int)transformMatrices.size());
    }

    void pop() {
        if(pushPoints.empty()) {
            cout << "Error: pop before push(pop ignored)" << endl;
        } else {
            //printMatrix(transformMatrices.top());
            int pushPoint = pushPoints.top();
            while(transformMatrices.size() != pushPoint) {
                transformMatrices.pop();
            }
            pushPoints.pop();
            //printMatrix(transformMatrices.top());
        }
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
        return new Triangle(transform(triangle->a, transformMatrix), transform(triangle->b, transformMatrix), transform(triangle->c, transformMatrix));
    }
};

int main() {
    Vector *eye, *look, *up;
    double fovY, aspectRatio, near, far;

    string sceneFileName = "scene.txt";
    string modelOutputFileName = "out_stage1.txt";
    string viewOutputFileName = "out_stage2.txt";
    string projectionOutputFileName = "out_stage3.txt";

    string line;
    ifstream sceneFile;
    ofstream outputFile;
    sceneFile.open(sceneFileName);

    getline(sceneFile, line);
    eye = new Vector(line);

    getline(sceneFile, line);
    look = new Vector(line);

    getline(sceneFile, line);
    up = new Vector(line);

    sceneFile >> fovY;
    sceneFile >> aspectRatio;
    sceneFile >> near;
    sceneFile >> far;

    GL gl;
    stack<double**> s;
    stack<int> pushPoints;
    vector<Triangle*> triangles;

    s.push(initializeMatrix());

    while(getline(sceneFile, line)) {
        if (line == "triangle") {
            Vector *a, *b, *c;
            getline(sceneFile, line);
            a = new Vector(line);

            getline(sceneFile, line);
            b = new Vector(line);

            getline(sceneFile, line);
            c = new Vector(line);

            //triangles.push_back(transformTriangle(new Triangle(a, b, c), s.top()));
            triangles.push_back(gl.transformTriangle(new Triangle(a, b, c), s.top()));
        }
        else if (line == "translate") {
            /*double** translateMatrix = initializeMatrix();
            sceneFile >> translateMatrix[0][3];
            sceneFile >> translateMatrix[1][3];
            sceneFile >> translateMatrix[2][3];*/

            double tx, ty, tz;
            sceneFile >> tx;
            sceneFile >> ty;
            sceneFile >> tz;

            s.push(gl.product(s.top(), gl.makeTranslateMatrix(tx, ty, tz)));
        }
        else if (line == "rotate") {
            double angle, x, y, z;

            sceneFile >> angle;
            sceneFile >> x;
            sceneFile >> y;
            sceneFile >> z;

            auto* axis = new Vector(x, y, z);
            axis->normalize();

            s.push(gl.product(s.top(), gl.makeRotationMatrix(axis, angle)));
        }
        else if (line == "scale") {
            double sx, sy, sz;
            sceneFile >> sx;
            sceneFile >> sy;
            sceneFile >> sz;

            s.push(product(s.top(), gl.makeScaleMatrix(sx, sy, sz)));
        }
        else if (line == "push") {
            pushPoints.push(s.size());

            //cout << "pushed with size " << s.size() << endl;
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
    sceneFile.close();

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
        outputFile << gl.transformTriangle(triangle, viewMatrix)->toString() << endl << endl;
    }
    outputFile.close();

    //task 3
    double** projectionMatrix = getProjectionMatrix(fovY, aspectRatio, near, far, viewMatrix);

    gl.printMatrix(projectionMatrix);

    vector<Triangle*> projectedTriangles;
    projectedTriangles.reserve(triangles.size());
    for (auto & triangle : triangles) {
        projectedTriangles.push_back(gl.transformTriangle(triangle, projectionMatrix));
    }

    outputFile.open(projectionOutputFileName);
    for (auto & triangle : projectedTriangles) {
        outputFile << triangle->toString() << endl << endl;
    }
    outputFile.close();



    return 0;
}
