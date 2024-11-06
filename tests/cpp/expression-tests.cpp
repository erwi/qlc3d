#include <catch.h>
#include <expression.h>
#include <geom/vec3.h>

TEST_CASE("Expression with no coordinates evaluates to same value everywhere") {
  CartesianExpression expr("3 + 3");
  REQUIRE(expr.evaluate(0, 0, 0) == 6);
  REQUIRE(expr.evaluate(Vec3(1, 2, 3)) == 6);
}

TEST_CASE("Expression with coordinates evaluates to different values") {
  CartesianExpression expr("X");
  REQUIRE(expr.evaluate(0, 0, 0) == 0);
  REQUIRE(expr.evaluate(Vec3(1, 2, 3)) == 1);
}

TEST_CASE("Expression with upper and lowercase x and Y work") {
  CartesianExpression expr1("x + Y");
  CartesianExpression expr2("X + y");

  REQUIRE(expr1.evaluate(1, 2, 0) == 3);
  REQUIRE(expr2.evaluate(1, 2, 0) == 3);
}

TEST_CASE("Throw exception with helpful information when expression is invalid") {
  try {
    CartesianExpression expr("1 + a"); // can't have a variable 'a' in the expression
  } catch (const ExpressionException &e) {
    REQUIRE(std::string(e.what()) == "Failed to compile expression: 1 + a. Error near location 5");
    REQUIRE(e.getLocation() == 5);
    return;
  }

  FAIL("Expected ExpressionException to be thrown");
}