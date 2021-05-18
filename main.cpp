#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

#define N 4
#define PI (2*acos(0.0))

class Point{
public:
    double x, y, z;

    Point() : x(0), y(0), z(0) {}

    Point(double x, double y, double z) : x(x), y(y), z(z) {}

    explicit Point(const string& line) : x(0), y(0), z(0) {
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

    string toString() const {
        return to_string(x) + " " + to_string(y) + " " + to_string(z);
    }
};

class Triangle{
public:
    Point *a, *b, *c;

    Triangle(Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}

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

Point* transform(double** matrix, const Point& p) {
    double resultMatrix[N];
    double pointMatrix[] = {p.x, p.y, p.z, 1};

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            resultMatrix[i] += matrix[i][j]*pointMatrix[j];
        }
    }

    return new Point(resultMatrix[0], resultMatrix[1], resultMatrix[2]);
}

Point rotate(Point* p, Point* axis, double angle) {
    Point crossProduct, result;
    double theta = angle * (PI/180);
    double dotProduct = p->x*axis->x + p->y*axis->y + p->z*axis->z;

    crossProduct.x = axis->y*p->z - axis->z*p->y;
    crossProduct.y = axis->z*p->x - axis->x*p->z;
    crossProduct.z = axis->x*p->y - axis->y*p->x;

    result.x = p->x * cos(theta) + crossProduct.x*sin(theta) + axis->x*dotProduct*(1-cos(theta));
    result.y = p->y * cos(theta) + crossProduct.y*sin(theta) + axis->y*dotProduct*(1-cos(theta));
    result.z = p->z * cos(theta) + crossProduct.z*sin(theta) + axis->z*dotProduct*(1-cos(theta));

    return result;
}

double** makeRotationMatrix(Point* axis, double angle) {
    Point c1 = rotate(new Point(1, 0, 0), axis, angle);
    Point c2 = rotate(new Point(0, 1, 0), axis, angle);
    Point c3 = rotate(new Point(0, 0, 1), axis, angle);

    double** result = initializeMatrix();
    result[0][0] = c1.x;
    result[1][0] = c1.y;
    result[2][0] = c1.z;

    result[0][1] = c2.x;
    result[1][1] = c2.y;
    result[2][1] = c2.z;

    result[0][2] = c3.x;
    result[1][2] = c3.y;
    result[2][2] = c3.z;

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
    //Point eye, look, up;
    double fovX = 1, fovY, aspectRatio, near, far;

    string data, line;
    string sceneFileName = "scene.txt";
    ifstream sceneFile;

    sceneFile.open(sceneFileName);

    getline(sceneFile, line);
    Point eye(line);

    getline(sceneFile, line);
    Point look(line);

    getline(sceneFile, line);
    Point up(line);

    sceneFile >> fovY;
    sceneFile >> aspectRatio;
    sceneFile >> near;
    sceneFile >> far;


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
            Point a(line);

            getline(sceneFile, line);
            Point b(line);

            getline(sceneFile, line);
            Point c(line);

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

            auto* axis = new Point(x, y, z);
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
