#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>

// NodeType enumeration
enum class NodeType {
    INNER,
    EDGE,
    DUMMY,
};

// Node structure
struct Node {
    NodeType type;
    double x;
    double y;
    Node(NodeType type, double x, double y) : type(type), x(x), y(y) {}
};

// Figure base class
class Figure {
public:
    virtual bool isEdgeNode(double x, double y, double eps) = 0;
    virtual bool isInnerNode(double x, double y, double eps) = 0;
    virtual bool isLeft(double x, double y, double eps) = 0;
    virtual bool isRight(double x, double y, double eps) = 0;
    virtual bool isBottom(double x, double y, double eps) = 0;
    virtual bool isTop(double x, double y, double eps) = 0;
};

// TShapedFigure class
class TShapedFigure : public Figure {
private:
    double _scale;
    double _ratio;

public:
    TShapedFigure(double scale, double ratio);
    bool isEdgeNode(double x, double y, double eps) override;
    bool isInnerNode(double x, double y, double eps) override;
    bool isLeft(double x, double y, double eps) override;
    bool isRight(double x, double y, double eps) override;
    bool isBottom(double x, double y, double eps) override;
    bool isTop(double x, double y, double eps) override;
};

TShapedFigure::TShapedFigure(double scale, double ratio) : _scale(scale), _ratio(ratio) {}

bool TShapedFigure::isEdgeNode(double x, double y, double eps) {
    double connectionHeight = _scale - _ratio * _scale;
    double widthOfT = _ratio * _scale;
    if (y < connectionHeight + eps && y > connectionHeight - eps) {
        return x < _scale / 2 - widthOfT / 2 + eps || x > _scale / 2 + widthOfT / 2 - eps;
    }
    if (y < connectionHeight + eps) {
        return ((y < 0.0 + eps && y > 0.0 - eps) && (x < _scale / 2 + widthOfT / 2 + eps && x > _scale / 2 - widthOfT / 2 - eps)) ||
            (x < _scale / 2 - widthOfT / 2 + eps && x > _scale / 2 - widthOfT / 2 - eps) ||
            (x < _scale / 2 + widthOfT / 2 + eps && x > _scale / 2 + widthOfT / 2 - eps);
    }
    return (x < 0.0 + eps && x > 0.0 - eps) || (x < _scale + eps && x > _scale - eps) || (y < _scale + eps && y > _scale - eps);
}

bool TShapedFigure::isInnerNode(double x, double y, double eps) {
    double connectionHeight = _scale - _ratio * _scale;
    double widthOfT = _ratio * _scale;
    if (y < connectionHeight + eps && y > connectionHeight - eps) {
        return x < _scale / 2 + widthOfT / 2 + eps && x > _scale / 2 - widthOfT / 2 - eps;
    }
    if (y < connectionHeight + eps) {
        return x < _scale / 2 + widthOfT / 2 + eps && x > _scale / 2 - widthOfT / 2 - eps && y > 0.0 - eps;
    }
    return y < _scale + eps && x > 0.0 + eps && x < _scale - eps;
}

bool TShapedFigure::isRight(double x, double y, double eps) {
    double connectionHeight = _scale - _ratio * _scale;
    double widthOfT = _ratio * _scale;
    if (y < connectionHeight + eps) {
        return (x < _scale / 2 + widthOfT / 2 + eps && x > _scale / 2 + widthOfT / 2 - eps);
    }
    return (x < _scale + eps && x > _scale - eps);
}
bool TShapedFigure::isLeft(double x, double y, double eps) {
    double connectionHeight = _scale - _ratio * _scale;
    double widthOfT = _ratio * _scale;
    if (y < connectionHeight + eps) {
        return (x < _scale / 2 - widthOfT / 2 + eps && x > _scale / 2 - widthOfT / 2 - eps);
    }
    return (x < 0.0 + eps && x > 0.0 - eps);
}

bool TShapedFigure::isBottom(double x, double y, double eps) {
    double connectionHeight = _scale - _ratio * _scale;
    double widthOfT = _ratio * _scale;
    return (y < 0.0 + eps && y > 0.0 - eps) && (x < _scale / 2 + widthOfT / 2 + eps && x > _scale / 2 - widthOfT / 2 - eps) ||
        (y < connectionHeight + eps && y > connectionHeight - eps);
}

bool TShapedFigure::isTop(double x, double y, double eps) {
    return (y < _scale + eps && y > _scale - eps);
}

// Square class
class Square : public Figure {
private:
    double _scale;

public:
    Square(double scale);
    bool isEdgeNode(double x, double y, double eps) override;
    bool isInnerNode(double x, double y, double eps) override;
    bool isLeft(double x, double y, double eps) override;
    bool isRight(double x, double y, double eps) override;
    bool isBottom(double x, double y, double eps) override;
    bool isTop(double x, double y, double eps) override;
};

Square::Square(double scale) : _scale(scale) {}

bool Square::isEdgeNode(double x, double y, double eps) {
    return (x < _scale + eps && x > _scale - eps) || (x < 0.0 + eps && x > 0.0 - eps) ||
        (y < _scale + eps && y > _scale - eps) || (y < 0.0 + eps && y > 0.0 - eps);
}

bool Square::isInnerNode(double x, double y, double eps) {
    return x < _scale + eps && y < _scale + eps && x > 0.0 - eps && y > 0.0 - eps;
}

bool Square::isLeft(double x, double y, double eps) {
    return (x < 0.0 + eps && x > 0.0 - eps);
}

bool Square::isRight(double x, double y, double eps) {
    return (x < _scale + eps && x > _scale - eps);
}

bool Square::isBottom(double x, double y, double eps) {
    return (y < 0.0 + eps && y > 0.0 - eps);
}

bool Square::isTop(double x, double y, double eps) {
    return (y < _scale + eps && y > _scale - eps);
}

// ThirdConditionsSide enumeration
enum class ThirdConditionsSide {
    LEFT,
    RIGHT,
    BOTTOM,
    TOP,
    NONE,
};

// RegularGrid class
class RegularGrid {
private:
    double _scale;
    double _stepX;
    double _stepY;
    Figure* _figure;
    std::function<double(double, double)> _u;
    std::function<double(double, double)> _f;
    ThirdConditionsSide _side;
    double _epsBase = 1e-5;
    double _beta = 1.0;
    double _lambda = 1.0;
    double _gamma = 0.0;

public:
    RegularGrid(double scale, double stepX, double stepY, Figure* figure,
        std::function<double(double, double)> f, std::function<double(double, double)> u, ThirdConditionsSide side);

    std::vector<Node> nodes();
    double leftDerivativeByX(double x, double y);
    double leftDerivativeByY(double x, double y);
    double rightDerivativeByX(double x, double y);
    double rightDerivativeByY(double x, double y);
};

RegularGrid::RegularGrid(double scale, double stepX, double stepY, Figure* figure,
    std::function<double(double, double)> f, std::function<double(double, double)> u, ThirdConditionsSide side)
    : _scale(scale), _stepX(stepX), _stepY(stepY), _figure(figure), _u(u), _f(f), _side(side) {
}

std::vector<Node> RegularGrid::nodes() {
    std::vector<Node> nodes;
    int numberOfNodesByX = static_cast<int>(_scale / _stepX);
    int numberOfNodesByY = static_cast<int>(_scale / _stepY);
    double eps = std::min(_stepX, _stepY) * _epsBase;

    for (int i = 0; i <= numberOfNodesByX; ++i) {
        for (int j = 0; j <= numberOfNodesByY; ++j) {
            double x = static_cast<double>(i) * _stepX;
            double y = static_cast<double>(j) * _stepY;
            if (_figure->isEdgeNode(x, y, eps)) {
                nodes.emplace_back(Node(NodeType::EDGE, x, y));
            }
            else if (_figure->isInnerNode(x, y, eps)) {
                nodes.emplace_back(Node(NodeType::INNER, x, y));
            }
            else {
                nodes.emplace_back(Node(NodeType::DUMMY, x, y));
            }
        }
    }
    return nodes;
}

double RegularGrid::leftDerivativeByX(double x, double y) {
    double h = 1e-9;
    return (_u(x, y) - _u(x - h, y)) / h;
}

double RegularGrid::leftDerivativeByY(double x, double y) {
    double h = 1e-9;
    return (_u(x, y) - _u(x, y - h)) / h;
}

double RegularGrid::rightDerivativeByX(double x, double y) {
    double h = 1e-9;
    return (_u(x + h, y) - _u(x, y)) / h;
}

double RegularGrid::rightDerivativeByY(double x, double y) {
    double h = 1e-9;
    return (_u(x, y + h) - _u(x, y)) / h;
}

// SystemOfEquations class
class SystemOfEquations {
protected:
    std::vector<double> _di, _au1, _au2, _al1, _al2;
    std::vector<double> _b;
    int _maxiter = 10000;
    double _w;
    double _eps = 1e-12;

public:
    SystemOfEquations(std::vector<double> di, std::vector<double> au1, std::vector<double> au2,
        std::vector<double> al1, std::vector<double> al2, std::vector<double> b);
    std::vector<double> solution();

protected:
    virtual std::vector<double> calculateXk(std::vector<double> x) = 0;
    double euclideanNorm(const std::vector<double>& x);
    double relativeDiscrepancy();
    std::vector<double> multiplyMatrixByVector(const std::vector<double>& x);
};

SystemOfEquations::SystemOfEquations(std::vector<double> di, std::vector<double> au1,
    std::vector<double> au2, std::vector<double> al1,
    std::vector<double> al2, std::vector<double> b)
    : _di(std::move(di)), _au1(std::move(au1)), _au2(std::move(au2)),
    _al1(std::move(al1)), _al2(std::move(al2)), _b(std::move(b)) {
}

std::vector<double> SystemOfEquations::solution() {
    int i = 0;
    std::vector<double> xk(_b.size(), 0.0);
    while (i < _maxiter && relativeDiscrepancy() >= _eps) {
        xk = calculateXk(xk);
        i++;
    }
    return xk;
}

double SystemOfEquations::euclideanNorm(const std::vector<double>& x) {
    double sum = 0.0;
    for (double val : x) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

double SystemOfEquations::relativeDiscrepancy() {
    int n = _b.size();
    std::vector<double> numerator(n);
    std::vector<double> multiplication = multiplyMatrixByVector(x);
    for (int i = 0; i < n; i++) {
        numerator[i] = _b[i] - multiplication[i];
    }
    return euclideanNorm(numerator) / euclideanNorm(_b);
}

std::vector<double> SystemOfEquations::multiplyMatrixByVector(const std::vector<double>& x) {
    std::vector<double> result(x.size(), 0.0);
    for (size_t i = 0; i < x.size(); ++i) {
        result[i] = _di[i] * x[i];
        if (i > 0) result[i] += _au1[i - 1] * x[i - 1];
        if (i < x.size() - 1) result[i] += _au2[i] * x[i + 1];
    }
    return result;
}

// Пример конкретной системы уравнений
class Type1System : public SystemOfEquations {
public:
    Type1System(const std::vector<double>& params)
        : SystemOfEquations(params[0], params[1], params[2], params[3], params[4], params[5]) {
    }

protected:
    std::vector<double> calculateXk(std::vector<double> x) override {
        // Пример итеративного метода, например, метод Гаусса-Зейделя
        for (size_t i = 0; i < x.size(); ++i) {
            double sum = _b[i];
            if (i > 0) sum -= _au1[i - 1] * x[i - 1];
            if (i < x.size() - 1) sum -= _au2[i] * x[i + 1];
            x[i] = sum / _di[i];
        }
        return x;
    }
};

class Type2System : public SystemOfEquations {
public:
    Type2System(const std::vector<double>& params)
        : SystemOfEquations(params[0], params[1], params[2], params[3], params[4], params[5]) {
    }

protected:
    std::vector<double> calculateXk(std::vector<double> x) override {
        // Пример другого итеративного метода
        for (size_t i = 0; i < x.size(); ++i) {
            double sum = _b[i];
            if (i > 0) sum -= _au1[i - 1] * x[i - 1];
            if (i < x.size() - 1) sum -= _au2[i] * x[i + 1];
            x[i] = sum / _di[i]; // Можно использовать другой алгоритм
        }
        return x;
    }
}

// Фабричный метод для создания систем уравнений
std::unique_ptr<SystemOfEquations> createSystem(const std::string& type, const std::vector<double>& params) {
    if (type == "type1") {
        return std::make_unique<Type1System>(params);
    }
    else if (type == "type2") {
        return std::make_unique<Type2System>(params);
    }
    throw std::invalid_argument("Unknown system type");
}

