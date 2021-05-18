#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <sstream>

using namespace std;

#define N 4

class Point{
public:
    double x, y, z;

    Point(double x, double y, double z) : x(x), y(y), z(z) {}

    Point(const string& line) {
        stringstream ss(line);
        ss >> x;
        ss >> y;
        ss >> z;
    }

    string toString() const {
        return to_string(x) + " " + to_string(y) + " " + to_string(z);
    }
};

class Triangle{
public:
    Point a, b, c;

    Triangle(const Point &a, const Point &b, const Point &c) : a(a), b(b), c(c) {}

    string toString() const {
        return a.toString() + "\n" + b.toString() + "\n" + c.toString();
    }
};

class Matrix{
public:
    int n;
    bool nop;
    double** matrix;

    explicit Matrix(int n, bool nop){
        this->n = n;
        this->nop = nop;
        matrix = new double*[n];

        for (int i = 0; i < n; ++i) {
            matrix[i] = new double[n];
            for (int j = 0; j < n; ++j) {
                matrix[i][j] = (i==j) ? 1: 0;
            }
        }
    }

    double** product(const Matrix& mat) const {
        double** result;
        result = new double*[n];

        for (int i = 0; i < n; ++i) {
            result[i] = new double[n];
            for (int j = 0; j < n; ++j) {
                result[i][j] = 0;
            }
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    result[i][j] += matrix[i][k]*mat.matrix[k][j];
                }
            }
        }

        return result;
    }

    Point* transform(const Point& p) const {
        double result[n];
        double pointMatrix[] = {p.x, p.y, p.z, 1};

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result[i] += matrix[i][j]*pointMatrix[j];
            }
        }

        return new Point(result[0], result[1], result[2]);
    }

    virtual ~Matrix() {
        for (int i = 0; i < n; ++i) {
            delete[] matrix[i];
        }
        delete[] matrix;
    }
};

double** product(double a[N][N], double b[N][N]) {
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


int main() {
    stack<Matrix*> s;
    vector<Triangle*> triangles;

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


    auto* top = new Matrix(N, false);
    s.push(top);

    while(getline(sceneFile, line)) {
        if (line == "triangle") {
            getline(sceneFile, line);
            Point a(line);

            getline(sceneFile, line);
            Point b(line);

            getline(sceneFile, line);
            Point c(line);

            triangles.push_back(new Triangle(*s.top()->transform(a), *s.top()->transform(a), *s.top()->transform(a)));
        }
        else if (line == "translate") {
            getline(sceneFile, line);
            Point translate(line);

            Matrix translateMatrix(N, false);
            translateMatrix.matrix[N-1][0] = translate.x;
            translateMatrix.matrix[N-1][1] = translate.y;
            translateMatrix.matrix[N-1][2] = translate.z;
        }
        else if (line == "rotate") {
            getline(sceneFile, line);
            //Point rotate(line);

        }
        else if (line == "scale") {
            getline(sceneFile, line);
            Point scale(line);

        }
        else if (line == "push") {

        }
        else if (line == "pop") {

        }
        else if (line == "end") {
            break;
        }
        else {
            cout << "something is wrong" << endl;
        }

    }
    sceneFile.close();

    for (auto & triangle : triangles) {
        cout << triangle->toString() << endl << endl;
    }

    return 0;
}
