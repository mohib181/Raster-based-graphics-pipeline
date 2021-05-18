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

    string toString() const {
        return to_string(x) + " " + to_string(y) + " " + to_string(z);
    }
};

class Triangle{
public:
    Vector *a, *b, *c;

    Triangle(Vector *a, Vector *b, Vector *c) : a(a), b(b), c(c) {}

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

Vector* transform(double** matrix, const Vector& p) {
    double resultMatrix[N];
    double pointMatrix[] = {p.x, p.y, p.z, 1};

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            resultMatrix[i] += matrix[i][j]*pointMatrix[j];
        }
    }

    return new Vector(resultMatrix[0], resultMatrix[1], resultMatrix[2]);
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

void printMatrix(double** matrix) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


int main() {
    //Vector eye, look, up;
    double fovX = 1, fovY, aspectRatio, near, far;

    string data, line;
    string sceneFileName = "scene.txt";
    ifstream sceneFile;

    sceneFile.open(sceneFileName);

    getline(sceneFile, line);
    Vector eye(line);

    getline(sceneFile, line);
    Vector look(line);

    getline(sceneFile, line);
    Vector up(line);

    sceneFile >> fovY;
    sceneFile >> aspectRatio;
    sceneFile >> near;
    sceneFile >> far;

    Vector l, r, u;
    l.x = look.x - eye.x;
    l.y = look.y - eye.y;
    l.z = look.z - eye.z;

    l.normalize();


    /*cout << eye.toString() << endl;
    cout << look.toString() << endl;
    cout << up.toString() << endl;

    cout << fovX << " " << fovY << endl;
    cout << aspectRatio << endl;
    cout << near << " " << far << endl;
    */


    //auto* top = new Matrix(N, false);
    stack<double**> s;
    stack<int> pushPoints;
    vector<Triangle*> triangles;
    cout << pushPoints.size() << endl;

    s.push(initializeMatrix());

    while(getline(sceneFile, line)) {
        if (line == "triangle") {
            getline(sceneFile, line);
            Vector a(line);

            getline(sceneFile, line);
            Vector b(line);

            getline(sceneFile, line);
            Vector c(line);

            //triangles.push_back(new Triangle(*s.top()->transform(a), *s.top()->transform(a), *s.top()->transform(a)));
            triangles.push_back(new Triangle(transform(s.top(), a), transform(s.top(), b), transform(s.top(), c)));
        }
        else if (line == "translate") {
            double** translateMatrix = initializeMatrix();
            sceneFile >> translateMatrix[0][3];
            sceneFile >> translateMatrix[1][3];
            sceneFile >> translateMatrix[2][3];

            cout << "translate" << endl;
            printMatrix(s.top());
            printMatrix(translateMatrix);
            s.push(product(s.top(), translateMatrix));
            printMatrix(s.top());
        }
        else if (line == "rotate") {
            double angle, x, y, z;

            sceneFile >> angle;
            sceneFile >> x;
            sceneFile >> y;
            sceneFile >> z;

            auto* axis = new Vector(x, y, z);
            axis->normalize();

            double** rotationMatrix = makeRotationMatrix(axis, angle);

            cout << "rotate" << endl;
            printMatrix(s.top());
            printMatrix(rotationMatrix);
            s.push(product(s.top(), rotationMatrix));
            printMatrix(s.top());
        }
        else if (line == "scale") {
            double** scaleMatrix = initializeMatrix();
            sceneFile >> scaleMatrix[0][0];
            sceneFile >> scaleMatrix[1][1];
            sceneFile >> scaleMatrix[2][2];

            cout << "scale" << endl;
            printMatrix(s.top());
            printMatrix(scaleMatrix);
            s.push(product(s.top(), scaleMatrix));
            printMatrix(s.top());

        }
        else if (line == "push") {
            pushPoints.push(s.size());
            cout << "pushed with size " << s.size() << endl;
            printMatrix(s.top());
        }
        else if (line == "pop") {
            if(pushPoints.empty()) {
                cout << "Error: pop before push(pop ignored)" << endl;
            } else {
                printMatrix(s.top());
                int pushPoint = pushPoints.top();
                while(s.size() != pushPoint) {
                    s.pop();
                }
                pushPoints.pop();
                printMatrix(s.top());
            }
        }
        else if (line == "end") {
            break;
        }

    }
    sceneFile.close();

    for (auto & triangle : triangles) {
        cout << triangle->toString() << endl << endl;
    }

    return 0;
}
