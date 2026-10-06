/**
 * @file expression.cpp
 * @brief Expression: parsing and evaluation of the derived-parameter expressions.
 */
#include "expression.h"

#include <cctype>
#include <cmath>
#include <cstdlib>

#include "logging.h"

namespace pytec3d
{
    namespace
    {
        using libMesh::Real;
        using Vec = libMesh::VectorValue<Real>;
        using Ten = libMesh::TensorValue<Real>;

        [[noreturn]] void fail(const std::string &where, const std::string &message)
        {
            logging::error_and_abort("Expression " + where + ": " + message);
        }

        [[noreturn]] void type_error(const std::string &where, const char *what, ParameterType a, ParameterType b)
        {
            fail(where, std::string(what) + " is not defined for " + type_name(a) + " and " + type_name(b) + ".");
        }

        void require(const std::string &where, const char *what, const ParameterValue &value, ParameterType type)
        {
            if (value.type != type)
                fail(where, std::string(what) + " needs a " + type_name(type) + ", got a " + type_name(value.type) + ".");
        }

        ParameterValue make_vector(const Vec &v)
        {
            ParameterValue out;
            out.type = ParameterType::Vector;
            out.vector = v;
            return out;
        }

        ParameterValue make_tensor(const Ten &t)
        {
            ParameterValue out;
            out.type = ParameterType::Tensor;
            out.tensor = t;
            return out;
        }

        ParameterValue add(const ParameterValue &a, const ParameterValue &b, Real sign, const std::string &where)
        {
            if (a.type != b.type)
                type_error(where, sign > 0 ? "addition" : "subtraction", a.type, b.type);
            ParameterValue out = a;
            out.scalar = a.scalar + sign * b.scalar;
            out.vector = Vec(a.vector + sign * b.vector);
            out.tensor = Ten(a.tensor + sign * b.tensor);
            return out;
        }

        ParameterValue multiply(const ParameterValue &a, const ParameterValue &b, const std::string &where)
        {
            if (a.type == ParameterType::Scalar)
                return scaled(b, a.scalar);
            if (b.type == ParameterType::Scalar)
                return scaled(a, b.scalar);
            if (a.type == ParameterType::Tensor && b.type == ParameterType::Vector)
                return make_vector(Vec(a.tensor * b.vector));
            if (a.type == ParameterType::Tensor && b.type == ParameterType::Tensor)
                return make_tensor(Ten(a.tensor * b.tensor));
            type_error(where, "multiplication", a.type, b.type);
        }

        // Applies f to a scalar, or to each entry of a vector.
        template <typename F>
        ParameterValue elementwise(const ParameterValue &value, F f, const char *what, const std::string &where)
        {
            if (value.type == ParameterType::Scalar)
                return scalar_value(f(value.scalar));
            if (value.type == ParameterType::Vector)
                return make_vector(Vec(f(value.vector(0)), f(value.vector(1)), f(value.vector(2))));
            fail(where, std::string(what) + " needs a scalar or a vector.");
        }

        // Inverse of a 3x3 tensor by cofactors; a diagonal tensor inverts entry by entry.
        // With positive_only a diagonal entry that is not positive and finite, or a full tensor with a
        // non-positive determinant, gives 0.
        Ten invert_tensor(const Ten &t, bool positive_only, const std::string &where)
        {
            const auto invert_entry = [&](Real x)
            {
                if (positive_only)
                    return (x > 0.0 && std::isfinite(x)) ? 1.0 / x : 0.0;
                return 1.0 / x;
            };

            bool diagonal = true;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = 0; j < 3; ++j)
                    if (i != j && t(i, j) != 0.0)
                        diagonal = false;
            if (diagonal)
                return Ten(invert_entry(t(0, 0)), 0.0, 0.0, 0.0, invert_entry(t(1, 1)), 0.0, 0.0, 0.0, invert_entry(t(2, 2)));

            const Real det = t.det();
            if (positive_only && !(det > 0.0 && std::isfinite(det)))
                return Ten();
            if (det == 0.0)
                fail(where, "inv of a singular tensor.");
            Ten out;
            out(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) / det;
            out(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) / det;
            out(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) / det;
            out(1, 0) = (t(1, 2) * t(2, 0) - t(1, 0) * t(2, 2)) / det;
            out(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) / det;
            out(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) / det;
            out(2, 0) = (t(1, 0) * t(2, 1) - t(1, 1) * t(2, 0)) / det;
            out(2, 1) = (t(0, 1) * t(2, 0) - t(0, 0) * t(2, 1)) / det;
            out(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) / det;
            return out;
        }
    } // namespace

    // ---- parsing ---------------------------------------------------------------------------------------------

    class Expression::Builder
    {
    public:
        Builder(Expression &expression, const std::string &text, const Resolver &resolve, unsigned int dimension, bool axisymmetric)
            : expression_(expression), text_(text), resolve_(resolve), dimension_(dimension), axisymmetric_(axisymmetric)
        {
        }

        void run()
        {
            expression_.root_ = parse_sum();
            skip_spaces();
            if (pos_ != text_.size())
                error("unexpected '" + text_.substr(pos_, 1) + "'");
        }

    private:
        [[noreturn]] void error(const std::string &message) const
        {
            fail(expression_.description_, message + " in '" + text_ + "'.");
        }

        void skip_spaces()
        {
            while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_])))
                ++pos_;
        }

        bool accept(char c)
        {
            skip_spaces();
            if (pos_ < text_.size() && text_[pos_] == c)
            {
                ++pos_;
                return true;
            }
            return false;
        }

        int add_node(const Node &node)
        {
            expression_.nodes_.push_back(node);
            return static_cast<int>(expression_.nodes_.size()) - 1;
        }

        int constant(const ParameterValue &value)
        {
            Node node;
            node.op = Op::Constant;
            node.constant = value;
            return add_node(node);
        }

        int binary(Op op, int a, int b)
        {
            Node node;
            node.op = op;
            node.a = a;
            node.b = b;
            return add_node(node);
        }

        int parse_sum()
        {
            int left = parse_product();
            for (;;)
            {
                if (accept('+'))
                    left = binary(Op::Add, left, parse_product());
                else if (accept('-'))
                    left = binary(Op::Subtract, left, parse_product());
                else
                    return left;
            }
        }

        int parse_product()
        {
            int left = parse_unary();
            for (;;)
            {
                if (accept('*'))
                    left = binary(Op::Multiply, left, parse_unary());
                else if (accept('/'))
                    left = binary(Op::Divide, left, parse_unary());
                else
                    return left;
            }
        }

        int parse_unary()
        {
            if (accept('-'))
            {
                Node node;
                node.op = Op::Negate;
                node.a = parse_unary();
                return add_node(node);
            }
            return parse_primary();
        }

        int parse_primary()
        {
            skip_spaces();
            if (pos_ >= text_.size())
                error("unexpected end");

            const char c = text_[pos_];
            if (c == '(')
            {
                ++pos_;
                const int inner = parse_sum();
                if (!accept(')'))
                    error("missing ')'");
                return inner;
            }
            if (std::isdigit(static_cast<unsigned char>(c)) || c == '.')
            {
                char *end = nullptr;
                const Real number = std::strtod(text_.c_str() + pos_, &end);
                if (end == text_.c_str() + pos_)
                    error("invalid number");
                pos_ = static_cast<std::size_t>(end - text_.c_str());
                return constant(scalar_value(number));
            }
            if (std::isalpha(static_cast<unsigned char>(c)) || c == '_')
                return parse_name();
            error(std::string("unexpected '") + c + "'");
        }

        int parse_name()
        {
            const std::size_t start = pos_;
            while (pos_ < text_.size() && (std::isalnum(static_cast<unsigned char>(text_[pos_])) || text_[pos_] == '_'))
                ++pos_;
            const std::string name = text_.substr(start, pos_ - start);

            // A parameter keeps its type suffix: name[][], name[], name().
            for (const char *suffix : {"[][]", "[]", "()"})
            {
                const std::string s(suffix);
                if (text_.compare(pos_, s.size(), s) == 0)
                {
                    pos_ += s.size();
                    return symbol(name + s);
                }
            }

            if (pos_ < text_.size() && text_[pos_] == '(')
                return parse_call(name);
            return parse_plain(name);
        }

        int symbol(const std::string &key)
        {
            const std::optional<unsigned int> index = resolve_(key);
            if (!index)
                error("unknown parameter '" + key + "'");
            Node node;
            node.op = Op::Symbol;
            node.symbol = *index;
            return add_node(node);
        }

        int parse_plain(const std::string &name)
        {
            if (name == "pi")
                return constant(scalar_value(libMesh::pi));
            if (name == "mu0")
                return constant(scalar_value(4.0 * libMesh::pi * 1e-7));
            if (name == "dimension")
                return constant(scalar_value(static_cast<Real>(dimension_)));
            if (name == "axisymmetric")
                return constant(scalar_value(axisymmetric_ ? 1.0 : 0.0));
            if (name == "ex")
                return constant(make_vector(Vec(1.0, 0.0, 0.0)));
            if (name == "ey")
                return constant(make_vector(Vec(0.0, 1.0, 0.0)));
            if (name == "ez")
                return constant(make_vector(Vec(0.0, 0.0, 1.0)));
            return symbol(name);
        }

        int parse_call(const std::string &name)
        {
            struct Entry
            {
                const char *name;
                Function function;
                unsigned int arity;
            };
            static const Entry functions[] = {
                {"abs", Function::Abs, 1},
                {"sqrt", Function::Sqrt, 1},
                {"cbrt", Function::Cbrt, 1},
                {"det", Function::Det, 1},
                {"diag", Function::Diag, 1},
                {"outer", Function::Outer, 2},
                {"dot", Function::Dot, 2},
                {"inv", Function::Inv, 1},
                {"invpos", Function::Invpos, 1},
                {"global", Function::Global, 1},
                {"select", Function::Select, 3},
            };

            const Entry *entry = nullptr;
            for (const Entry &candidate : functions)
                if (name == candidate.name)
                    entry = &candidate;
            if (!entry)
                error("unknown function '" + name + "'");

            accept('(');
            int args[3] = {-1, -1, -1};
            for (unsigned int i = 0; i < entry->arity; ++i)
            {
                if (i > 0 && !accept(','))
                    error("'" + name + "' needs " + std::to_string(entry->arity) + " arguments");
                args[i] = parse_sum();
            }
            if (!accept(')'))
                error("missing ')' after '" + name + "' arguments");

            Node node;
            node.op = Op::Call;
            node.function = entry->function;
            node.a = args[0];
            node.b = args[1];
            node.c = args[2];
            return add_node(node);
        }

        Expression &expression_;
        const std::string &text_;
        const Resolver &resolve_;
        unsigned int dimension_;
        bool axisymmetric_;
        std::size_t pos_ = 0;
    };

    Expression Expression::parse(const std::string &text,
                                 const std::string &context,
                                 const Resolver &resolve,
                                 unsigned int dimension,
                                 bool axisymmetric)
    {
        Expression expression;
        expression.description_ = context;
        Builder(expression, text, resolve, dimension, axisymmetric).run();
        return expression;
    }

    // ---- evaluation ------------------------------------------------------------------------------------------

    ParameterValue Expression::evaluate(ExpressionContext &context) const
    {
        return evaluate_node(root_, context);
    }

    ParameterValue Expression::evaluate_node(int index, ExpressionContext &context) const
    {
        const Node &node = nodes_[index];
        const std::string &where = description_;

        switch (node.op)
        {
        case Op::Constant:
            return node.constant;
        case Op::Symbol:
            return context.symbol(node.symbol);
        case Op::Negate:
            return scaled(evaluate_node(node.a, context), -1.0);
        case Op::Add:
            return add(evaluate_node(node.a, context), evaluate_node(node.b, context), 1.0, where);
        case Op::Subtract:
            return add(evaluate_node(node.a, context), evaluate_node(node.b, context), -1.0, where);
        case Op::Multiply:
            return multiply(evaluate_node(node.a, context), evaluate_node(node.b, context), where);
        case Op::Divide:
        {
            const ParameterValue numerator = evaluate_node(node.a, context);
            const ParameterValue denominator = evaluate_node(node.b, context);
            require(where, "division", denominator, ParameterType::Scalar);
            return scaled(numerator, 1.0 / denominator.scalar);
        }
        case Op::Call:
            break;
        }

        // select() evaluates the branch it takes only: the other one may use parameters the element does not define.
        if (node.function == Function::Select)
        {
            const ParameterValue condition = evaluate_node(node.a, context);
            require(where, "select condition", condition, ParameterType::Scalar);
            return evaluate_node(condition.scalar != 0.0 ? node.b : node.c, context);
        }

        const ParameterValue a = evaluate_node(node.a, context);
        const ParameterValue b = node.b >= 0 ? evaluate_node(node.b, context) : ParameterValue();
        const ParameterValue c = node.c >= 0 ? evaluate_node(node.c, context) : ParameterValue();

        switch (node.function)
        {
        case Function::Abs:
            require(where, "abs", a, ParameterType::Scalar);
            return scalar_value(std::abs(a.scalar));
        case Function::Sqrt:
            require(where, "sqrt", a, ParameterType::Scalar);
            return scalar_value(std::sqrt(a.scalar));
        case Function::Cbrt:
            require(where, "cbrt", a, ParameterType::Scalar);
            return scalar_value(std::cbrt(a.scalar));
        case Function::Det:
            require(where, "det", a, ParameterType::Tensor);
            return scalar_value(a.tensor.det());
        case Function::Diag:
            require(where, "diag", a, ParameterType::Vector);
            return make_tensor(Ten(a.vector(0), 0.0, 0.0, 0.0, a.vector(1), 0.0, 0.0, 0.0, a.vector(2)));
        case Function::Outer:
        {
            require(where, "outer", a, ParameterType::Vector);
            require(where, "outer", b, ParameterType::Vector);
            Ten t;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = 0; j < 3; ++j)
                    t(i, j) = a.vector(i) * b.vector(j);
            return make_tensor(t);
        }
        case Function::Dot:
        {
            require(where, "dot", a, ParameterType::Vector);
            require(where, "dot", b, ParameterType::Vector);
            return scalar_value(a.vector(0) * b.vector(0) + a.vector(1) * b.vector(1) + a.vector(2) * b.vector(2));
        }
        case Function::Inv:
            if (a.type == ParameterType::Tensor)
                return make_tensor(invert_tensor(a.tensor, false, where));
            return elementwise(a, [](Real x)
                               { return 1.0 / x; }, "inv", where);
        case Function::Invpos:
            if (a.type == ParameterType::Tensor)
                return make_tensor(invert_tensor(a.tensor, true, where));
            return elementwise(a, [](Real x)
                               { return (x > 0.0 && std::isfinite(x)) ? 1.0 / x : 0.0; }, "invpos", where);
        case Function::Global:
        {
            const Ten &R = context.rotation();
            if (a.type == ParameterType::Vector)
                return make_vector(Vec(R * a.vector));
            if (a.type == ParameterType::Tensor)
                return make_tensor(Ten(R * a.tensor * R.transpose()));
            fail(where, "global needs a vector or a tensor.");
        }
        case Function::Select:
            break; // Handled above.
        }
        fail(where, "unknown function.");
    }
} // namespace pytec3d
