/**
 * @file expression.cpp
 * @brief Expression: recursive-descent parser to a stack machine, and its evaluation.
 */
#include "expression.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <utility>

namespace pytec3d
{
    namespace
    {
        enum UnaryOp
        {
            Negate,
            Not
        };
        enum BinaryOp
        {
            Add,
            Subtract,
            Multiply,
            Divide,
            Power,
            Less,
            LessEqual,
            Greater,
            GreaterEqual,
            Equal,
            NotEqual,
            And,
            Or
        };
        enum FunctionId
        {
            Sin,
            Cos,
            Tan,
            Asin,
            Acos,
            Atan,
            Atan2,
            Sinh,
            Cosh,
            Tanh,
            Exp,
            Log,
            Log10,
            Sqrt,
            Abs,
            Sign,
            Floor,
            Ceil,
            Pow,
            Min,
            Max,
            Hypot
        };

        /** @brief Name -> (id, number of arguments; -1: two or more). */
        const std::map<std::string, std::pair<int, int>> &functions()
        {
            static const std::map<std::string, std::pair<int, int>> table = {
                {"sin", {Sin, 1}},     {"cos", {Cos, 1}},     {"tan", {Tan, 1}},     {"asin", {Asin, 1}},
                {"acos", {Acos, 1}},   {"atan", {Atan, 1}},   {"atan2", {Atan2, 2}}, {"sinh", {Sinh, 1}},
                {"cosh", {Cosh, 1}},   {"tanh", {Tanh, 1}},   {"exp", {Exp, 1}},     {"log", {Log, 1}},
                {"log10", {Log10, 1}}, {"sqrt", {Sqrt, 1}},   {"abs", {Abs, 1}},     {"sign", {Sign, 1}},
                {"floor", {Floor, 1}}, {"ceil", {Ceil, 1}},   {"pow", {Pow, 2}},     {"min", {Min, -1}},
                {"max", {Max, -1}},    {"hypot", {Hypot, 2}}};
            return table;
        }

        const std::map<std::string, double> &constants()
        {
            static const std::map<std::string, double> table = {
                {"pi", 3.14159265358979323846}, {"e", 2.71828182845904523536}, {"mu0", 4e-7 * 3.14159265358979323846}};
            return table;
        }
    } // namespace

    /** @brief Builds the code of an Expression (one parse). */
    class ExpressionParser
    {
    public:
        ExpressionParser(Expression &out) : out_(out), s_(out.text_) {}

        void parse()
        {
            conditional();
            skip();
            if (pos_ != s_.size())
                fail("unexpected '" + std::string(1, s_[pos_]) + "'");
            out_.stack_size_ = max_depth_;
        }

    private:
        [[noreturn]] void fail(const std::string &what) const
        {
            throw std::invalid_argument("expression '" + s_ + "', at character " + std::to_string(pos_ + 1) + ": " + what +
                                        ".");
        }
        void skip()
        {
            while (pos_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[pos_])))
                ++pos_;
        }
        bool accept(const char *token)
        {
            skip();
            const std::size_t n = std::char_traits<char>::length(token);
            if (s_.compare(pos_, n, token) == 0)
            {
                pos_ += n;
                return true;
            }
            return false;
        }
        void emit(Expression::Op op, int pushed)
        {
            out_.code_.push_back(op);
            depth_ += pushed;
            max_depth_ = std::max(max_depth_, static_cast<std::size_t>(std::max(depth_, 0)));
        }
        void binary(int code) { emit({Expression::Op::Binary, code, 0, 0.0}, -1); }

        void conditional()
        {
            logical_or();
            if (accept("?"))
            {
                conditional();
                if (!accept(":"))
                    fail("':' expected");
                conditional();
                emit({Expression::Op::Select, 0, 0, 0.0}, -2);
            }
        }
        void logical_or()
        {
            logical_and();
            while (accept("||"))
            {
                logical_and();
                binary(Or);
            }
        }
        void logical_and()
        {
            comparison();
            while (accept("&&"))
            {
                comparison();
                binary(And);
            }
        }
        void comparison()
        {
            sum();
            static const std::pair<const char *, int> ops[] = {{"<=", LessEqual}, {">=", GreaterEqual}, {"==", Equal},
                                                               {"!=", NotEqual},  {"<", Less},          {">", Greater}};
            for (const auto &[token, code] : ops)
                if (accept(token))
                {
                    sum();
                    binary(code);
                    return;
                }
        }
        void sum()
        {
            product();
            while (true)
            {
                if (accept("+"))
                {
                    product();
                    binary(Add);
                }
                else if (accept("-"))
                {
                    product();
                    binary(Subtract);
                }
                else
                    return;
            }
        }
        void product()
        {
            unary();
            while (true)
            {
                skip();
                if (s_.compare(pos_, 2, "**") == 0)
                    return; // a power, handled by power()
                if (accept("*"))
                {
                    unary();
                    binary(Multiply);
                }
                else if (accept("/"))
                {
                    unary();
                    binary(Divide);
                }
                else
                    return;
            }
        }
        void unary()
        {
            if (accept("-"))
            {
                unary();
                emit({Expression::Op::Unary, Negate, 0, 0.0}, 0);
            }
            else if (accept("+"))
                unary();
            else if (accept("!"))
            {
                unary();
                emit({Expression::Op::Unary, Not, 0, 0.0}, 0);
            }
            else
                power();
        }
        void power()
        {
            primary();
            if (accept("^") || accept("**"))
            {
                unary(); // right associative: a ^ b ^ c = a ^ (b ^ c), a ^ -b allowed
                binary(Power);
            }
        }
        void primary()
        {
            skip();
            if (pos_ >= s_.size())
                fail("unexpected end");
            const char c = s_[pos_];
            if (std::isdigit(static_cast<unsigned char>(c)) || c == '.')
            {
                const char *begin = s_.c_str() + pos_;
                char *end = nullptr;
                const double value = std::strtod(begin, &end);
                if (end == begin)
                    fail("a number expected");
                pos_ += static_cast<std::size_t>(end - begin);
                emit({Expression::Op::Constant, 0, 0, value}, 1);
                return;
            }
            if (std::isalpha(static_cast<unsigned char>(c)) || c == '_')
            {
                const std::size_t start = pos_;
                while (pos_ < s_.size() && (std::isalnum(static_cast<unsigned char>(s_[pos_])) || s_[pos_] == '_'))
                    ++pos_;
                const std::string name = s_.substr(start, pos_ - start);
                if (accept("("))
                {
                    const auto it = functions().find(name);
                    if (it == functions().end())
                    {
                        pos_ = start;
                        fail("unknown function '" + name + "'");
                    }
                    int n = 0;
                    do
                    {
                        conditional();
                        ++n;
                    } while (accept(","));
                    if (!accept(")"))
                        fail("')' expected");
                    const int arity = it->second.second;
                    if ((arity > 0 && n != arity) || (arity < 0 && n < 2))
                    {
                        pos_ = start;
                        fail("'" + name + "' takes " + (arity > 0 ? std::to_string(arity) : std::string("two or more")) +
                             " argument(s), got " + std::to_string(n));
                    }
                    emit({Expression::Op::Function, it->second.first, n, 0.0}, 1 - n);
                    return;
                }
                const auto &vars = out_.variables_;
                const auto v = std::find(vars.begin(), vars.end(), name);
                if (v != vars.end())
                {
                    emit({Expression::Op::Variable, static_cast<int>(v - vars.begin()), 0, 0.0}, 1);
                    return;
                }
                const auto k = constants().find(name);
                if (k != constants().end())
                {
                    emit({Expression::Op::Constant, 0, 0, k->second}, 1);
                    return;
                }
                pos_ = start;
                std::string known;
                for (const auto &var : vars)
                    known += (known.empty() ? "" : ", ") + var;
                fail("unknown name '" + name + "' (the variables: " + (known.empty() ? "none" : known) +
                     "; the constants pi, e, mu0)");
            }
            if (accept("("))
            {
                conditional();
                if (!accept(")"))
                    fail("')' expected");
                return;
            }
            fail("unexpected '" + std::string(1, c) + "'");
        }

        Expression &out_;
        const std::string &s_;
        std::size_t pos_ = 0;
        int depth_ = 0;
        std::size_t max_depth_ = 0;
    };

    Expression::Expression(const std::string &text, std::vector<std::string> variables)
        : text_(text), variables_(std::move(variables))
    {
        ExpressionParser(*this).parse();
    }

    double Expression::evaluate(const double *values) const
    {
        // A small fixed stack for the usual expressions, a vector beyond.
        double small[32];
        std::vector<double> large;
        double *stack = small;
        if (stack_size_ > 32)
        {
            large.resize(stack_size_);
            stack = large.data();
        }
        std::size_t top = 0;
        for (const Op &op : code_)
        {
            switch (op.kind)
            {
            case Op::Constant: stack[top++] = op.value; break;
            case Op::Variable: stack[top++] = values[op.code]; break;
            case Op::Unary:
            {
                double &a = stack[top - 1];
                a = op.code == Negate ? -a : (a == 0.0 ? 1.0 : 0.0);
                break;
            }
            case Op::Binary:
            {
                const double b = stack[--top];
                double &a = stack[top - 1];
                switch (op.code)
                {
                case Add: a += b; break;
                case Subtract: a -= b; break;
                case Multiply: a *= b; break;
                case Divide: a /= b; break;
                case Power: a = std::pow(a, b); break;
                case Less: a = a < b; break;
                case LessEqual: a = a <= b; break;
                case Greater: a = a > b; break;
                case GreaterEqual: a = a >= b; break;
                case Equal: a = a == b; break;
                case NotEqual: a = a != b; break;
                case And: a = (a != 0.0) && (b != 0.0); break;
                case Or: a = (a != 0.0) || (b != 0.0); break;
                default: break;
                }
                break;
            }
            case Op::Select:
            {
                const double b = stack[--top];
                const double a = stack[--top];
                double &c = stack[top - 1];
                c = c != 0.0 ? a : b;
                break;
            }
            case Op::Function:
            {
                double *args = stack + top - op.arity;
                double r = 0.0;
                switch (op.code)
                {
                case Sin: r = std::sin(args[0]); break;
                case Cos: r = std::cos(args[0]); break;
                case Tan: r = std::tan(args[0]); break;
                case Asin: r = std::asin(args[0]); break;
                case Acos: r = std::acos(args[0]); break;
                case Atan: r = std::atan(args[0]); break;
                case Atan2: r = std::atan2(args[0], args[1]); break;
                case Sinh: r = std::sinh(args[0]); break;
                case Cosh: r = std::cosh(args[0]); break;
                case Tanh: r = std::tanh(args[0]); break;
                case Exp: r = std::exp(args[0]); break;
                case Log: r = std::log(args[0]); break;
                case Log10: r = std::log10(args[0]); break;
                case Sqrt: r = std::sqrt(args[0]); break;
                case Abs: r = std::abs(args[0]); break;
                case Sign: r = (args[0] > 0.0) - (args[0] < 0.0); break;
                case Floor: r = std::floor(args[0]); break;
                case Ceil: r = std::ceil(args[0]); break;
                case Pow: r = std::pow(args[0], args[1]); break;
                case Hypot: r = std::hypot(args[0], args[1]); break;
                case Min: r = *std::min_element(args, args + op.arity); break;
                case Max: r = *std::max_element(args, args + op.arity); break;
                default: break;
                }
                top -= static_cast<std::size_t>(op.arity);
                stack[top++] = r;
                break;
            }
            }
        }
        return top ? stack[top - 1] : 0.0;
    }
} // namespace pytec3d
