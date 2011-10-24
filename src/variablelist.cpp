#include <ctype.h>
#include "variablelist.h"

/*
 * Returns true if the given name already exists in the variable list
 */
bool Variablelist::exist(const char* name)
{
    return (get_id(name) != -1);
}


/*
 * Add a name and value to the variable list
 */
bool Variablelist::add(const char* name, double value)
{
    VAR new_var;
    strncpy(new_var.name, name, 30);
    new_var.value = value;

    int id = get_id(name);
    if (id == -1)
    {
        // variable does not yet exist
        var.push_back(new_var);
    }
    else
    {
        // variable already exists. overwrite it
        var[id] = new_var;
    }
    return true;
}

/*
 * Delete given variablename from the variable list
 */
bool Variablelist::del(const char* name)
{
    int id = get_id(name);
    if (id != -1)
    {
        var[id] = var[var.size()-1]; // move last item to deleted item
        var.pop_back();              // remove last item
        return true;
    }
    return false; 
}

/*
 * Get value of variable with given name
 */
bool Variablelist::get_value(const char* name, double* value)
{
    int id = get_id(name);
    if (id != -1)
    {
        *value = var[id].value;
        return true;
    }
    return false;
}


/*
 * Get value of variable with given id
 */
bool Variablelist::get_value(const int id, double* value)
{
    if (id >=0 && id < var.size())
    {
        *value = var[id].value;
        return true;
    }
    return false;
}


bool Variablelist::set_value(const char* name, const double value)
{
    return add(name, value);
}

/*
 * Returns the id of the given name in the variable list. Returns -1 if name
 * is not present in the list. Name is case insensitive
 */
int Variablelist::get_id(const char* name)
{
    // first make the name uppercase
    char nameU[NAME_LEN_MAX+1];
    char varU[NAME_LEN_MAX+1];
    toupper(nameU, name);
    
    for (int i = 0; i < var.size(); i++)
    {
        toupper(varU, var[i].name);
        if (strcmp(nameU, varU) == 0)
        {
            return i;
        }
    }
    return -1;
}


/*
 * str is copied to upper and made uppercase
 * upper is the returned string
 * str should be null-terminated
 */
void toupper(char upper[], const char str[])
{
    int i = -1;
    do
    {
        i++;
        upper[i] = toupper(str[i]);
    }
    while (str[i] != '\0');

}
