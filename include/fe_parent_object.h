#ifndef FEM_PARENT_OBJECT_H
#define FEM_PARENT_OBJECT_H

template <typename T1>
class FEParentObjectBase
{
public:
    FEParentObjectBase(T1 *parent = nullptr)
        : parent_object(parent) {}
    // Virtual destructor for proper cleanup of derived classes
    virtual ~FEParentObjectBase() = default;

    // Getter for the parent object
    T1 *get_parent_object() const
    {
        return parent_object;
    }
    // Setter for the parent object
    void set_parent_object(T1 *parent)
    {
        parent_object = parent;
    }

protected:
    // Attributes
    T1 *parent_object; // Stores the parent
};

#endif // FEM_PARENT_OBJECT_H