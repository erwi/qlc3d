#include <expression.h>
#include <geom/vec3.h>
#include <util/exception.h>
#include <util/stringutil.h>
#include <iostream>

CartesianExpression::CartesianExpression(const std::string &expression):
  expression_(StringUtil::toLowerCase(expression)),
  compiled_expression_(nullptr),
  x_(0), y_(0), z_(0)
{
  initialise();
}

CartesianExpression::CartesianExpression(const CartesianExpression &other):
  expression_(other.expression_),
  compiled_expression_(nullptr), // don't copy the compiled expression as it points to private variables x,y,z of the source object. It will be recompiled when needed.
  x_(other.x_), y_(other.y_), z_(other.z_)
{
  initialise();
}

void CartesianExpression::initialise() {
  int error = 0;
  if (compiled_expression_ == nullptr) {
    te_variable cartesianVariables[] = {{"x", &x_}, {"y", &y_}, {"z", &z_}};
    compiled_expression_ = te_compile(expression_.c_str(), cartesianVariables, 3, &error);
  } else {
    throw ExpressionException("Expression " + expression_ + " is already initialised", -1);
  }

  if (compiled_expression_ == nullptr) {
    throw ExpressionException("Failed to compile expression: " + expression_+ ". Error near location " + std::to_string(error), error);
  }
}

double CartesianExpression::evaluate(double x, double y, double z) const {
  x_ = x;
  y_ = y;
  z_ = z;
  return te_eval(compiled_expression_);
}

double CartesianExpression::evaluate(const Vec3& p) const {
  return evaluate(p.x(), p.y(), p.z());
}
