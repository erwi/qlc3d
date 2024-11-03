#include <expression.h>
#include <geom/vec3.h>
#include <util/exception.h>
#include <util/stringutil.h>

CartesianExpression::CartesianExpression(const std::string &expression):
  expression_(StringUtil::toLowerCase(expression)),
  compiled_expression_(nullptr)
{

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

double CartesianExpression::evaluate(double x, double y, double z) {
  if (compiled_expression_ == nullptr) {
    initialise();
  }

  x_ = x;
  y_ = y;
  z_ = z;
  return te_eval(compiled_expression_);
}
double CartesianExpression::evaluate(const Vec3& p) {
  return evaluate(p.x(), p.y(), p.z());
}

