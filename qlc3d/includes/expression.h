#ifndef PROJECT_QLC3D_EXPRESSION_H
#define PROJECT_QLC3D_EXPRESSION_H
#include <string>
#include <stdexcept>
#include <tinyexpr.h>

class Vec3;

class CartesianExpression {
  std::string expression_;
  te_expr *compiled_expression_;
  mutable double x_{}, y_{}, z_{};
  /* Initialise and compile the expression. Throws ExpressionException if the expression is invalid.*/
  void initialise();
public:
  explicit CartesianExpression(const std::string &expression);
  CartesianExpression(const CartesianExpression &other);

  [[nodiscard]] double evaluate(double x, double y, double z) const ; // TODO: try to make const, currently may compile and modifies x_, y_, z_
  [[nodiscard]] double evaluate(const Vec3 &p) const ;
  [[nodiscard]] const std::string &getExpression() const { return expression_; }
};

class ExpressionException : public std::runtime_error {
  int location_;
public:
  explicit ExpressionException(const std::string& message, int location) : std::runtime_error(message), location_(location) {}
  [[nodiscard]] int getLocation() const { return location_; }
};


#endif //PROJECT_QLC3D_EXPRESSION_H
