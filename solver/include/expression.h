/**
 * @file expression.h
 * @brief Expression: a real expression of named variables, parsed once from text and evaluated many times (the
 *        analytic curves of the nonlinear laws, model_data "nonlinear_laws").
 *
 * Grammar (the usual precedence; ^ is right associative and binds tighter than a unary minus, so -x^2 = -(x^2)):
 *
 *   expression  := conditional
 *   conditional := or ( "?" expression ":" expression )?
 *   or          := and ( "||" and )*
 *   and         := comparison ( "&&" comparison )*
 *   comparison  := sum ( ("<" | "<=" | ">" | ">=" | "==" | "!=") sum )?
 *   sum         := product ( ("+" | "-") product )*
 *   product     := unary ( ("*" | "/") unary )*
 *   unary       := ("-" | "+" | "!") unary | power
 *   power       := primary ( ("^" | "**") unary )?
 *   primary     := number | name | name "(" expression ("," expression)* ")" | "(" expression ")"
 *
 * Functions: sin cos tan asin acos atan atan2 sinh cosh tanh exp log log10 sqrt abs sign floor ceil pow min max
 * hypot. Constants: pi, e, mu0 (4 pi 1e-7). Comparisons and logical operators give 1 or 0. The variables are the
 * names given at parsing, in that order: evaluate(values) takes their values. A name that is neither a variable, a
 * constant nor a function is an error at parsing (with its position).
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

namespace pytec3d
{
    class Expression
    {
    public:
        Expression() = default;
        /** @brief Parses @p text with the variables @p variables; throws std::invalid_argument with the position of an
         * error. */
        Expression(const std::string &text, std::vector<std::string> variables);

        /** @brief The value at the variables' @p values (as many as variables, in their order). */
        double evaluate(const double *values) const;
        double operator()(double x) const { return evaluate(&x); }

        const std::string &text() const { return text_; }
        const std::vector<std::string> &variables() const { return variables_; }
        bool empty() const { return code_.empty(); }

    private:
        /** @brief One instruction of the stack machine the expression is compiled to. */
        struct Op
        {
            enum Kind : unsigned char
            {
                Constant,
                Variable,
                Unary,
                Binary,
                Function,
                Select, ///< c ? a : b (the three values on the stack)
            } kind;
            int code = 0;  ///< Variable index, operator or function id.
            int arity = 0; ///< Function arguments.
            double value = 0.0;
        };
        friend class ExpressionParser;

        std::string text_;
        std::vector<std::string> variables_;
        std::vector<Op> code_;
        std::size_t stack_size_ = 0;
    };
} // namespace pytec3d
