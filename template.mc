/**
 * @file template.mc
 * @author Alejandro Marquez
 * @date May 18, 2017
 * @brief Breve explicacion
 *
 * Here typically goes a more extensive explanation of what the header
 * defines. Doxygens tags are words preceeded by either a backslash @\
 * or by an at symbol @@.
 * @see http://www.stack.nl/~dimitri/doxygen/docblocks.html
 */
/**
 * @brief Ejemplo showing how to document a function with Doxygen.
 *
 * Description of what the function does. This part may refer to the parameters
 * of the function, like @p param1 or @p param2. A word of code can also be
 * inserted like @c this which is equivalent to <tt>this</tt> and can be useful
 * to say that the function returns a @c void or an @c int. If you want to have
 * more than one word in typewriter font, then just use @<tt@>.
 * We can also include text verbatim,
 * @verbatim like this@endverbatim
 * Sometimes it is also convenient to include an example of usage:
 * @code
 * systdef(eq,matrix([x[1](t)*u(t-1)],[u(t-2)]))$
 * eq_F;
 * eq_dF;
 * @endcode
 * Or,
 * @code{.py}
 * pyval = python_func(arg1, arg2)
 * print pyval
 * @endcode
 * when the language is not the one used in the current source file (but
 * <b>be careful</b> as this may be supported only by recent versions
 * of Doxygen). By the way, <b>this is how you write bold text</b> or,
 * if it is just one word, then you can just do @b this.
 * @param _name Descripción chida del first parameter of the function.
 * @param param2 The second one, which follows @p param1.
 * @return   Explicación de lo que regresa la función
 *  \f[
    |I_2|=\left| \int_{0}^T \psi(t)
             \left\{
                u(a,t)-
                \int_{\gamma(t)}^a
                \frac{d\theta}{k(\theta,t)}
                \int_{a}^\theta c(\xi)u_t(\xi,t)\,d\xi
             \right\} dt
          \right|
  \f]
 * @see Box_The_Last_One
 * @see http://website/
 * @note Something to note.
 * @warning Warning.
 */

systdef(
matrix _name,  /**< KKSome documentation for the member  symbol */
int _eq    /**< KKSome documentation for the member matrix */
)
{
  code...
};
