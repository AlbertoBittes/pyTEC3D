/**
 * @file fe_parent_object.h
 * @brief Non-owning parent-pointer mixin used throughout the solver hierarchy.
 *
 * Provides a lightweight link from a child object (assembler, builder, constraint
 * handler, etc.) back to its owning parent without creating hard include-cycle
 * dependencies.  The pointer is non-owning; lifetime management remains the
 * responsibility of the caller.
 */
#ifndef FEM_PARENT_OBJECT_H
#define FEM_PARENT_OBJECT_H

/**
 * @brief Mixin that stores a raw, non-owning pointer to a parent object.
 *
 * @tparam T1  Type of the parent (e.g. pyTEC3DEquationSystems, pyTEC3DLinearSystem).
 */
template <typename T1>
class FEParentObjectBase
{
public:
    explicit FEParentObjectBase(T1 *parent = nullptr)
        : parent_object(parent) {}

    virtual ~FEParentObjectBase() = default;

    /** @brief Returns the non-owning pointer to the parent object. */
    T1 *get_parent_object() const
    {
        return parent_object;
    }

    /** @brief Sets the non-owning pointer to the parent object. */
    void set_parent_object(T1 *parent)
    {
        parent_object = parent;
    }

protected:
    T1 *parent_object; ///< Raw non-owning pointer to the owning parent.
};

#endif // FEM_PARENT_OBJECT_H