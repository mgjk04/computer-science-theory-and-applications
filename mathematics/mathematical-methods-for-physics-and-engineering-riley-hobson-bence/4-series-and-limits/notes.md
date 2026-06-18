# Series and Limits

## Selected Concepts

### Definition of asymptotic approximation

For 2 functions $f(x)$ and $g(x)$ where $g(x) > 0$,
$$\lim_{x \to c}{\frac{f(x)}{g(x)}}=l,\ l \neq 0 \implies f(x) \approx  lg(x) $$
where $c$ is the value $x$ approaches.

### Taylor's theorem

For a function $f(x)$, assuming that for a given range of x, $f(x)$ is a continuous, single-valued function with derivatives $f'(x), f''(x),\ ...\ , f^{n}{(x)}.$
$$\int_{a}^{a + h}{f'(x)}\ dx = f(a + h) - f(a)$$
$$f(a + h) = f(a) + \int_{a}^{a + h}{f'(x)}\ dx$$
By substituting $f'(x) = f'(a)$, we may obtain the first approximation:
$$f(a + h) \approx f(a) + hf'(a)$$
In terms of $x$ and $a$ by substituting $x = a + h$,
$$f(x) \approx f(a) + (x - a)f'(a)$$
Similarly,
$$f'(x) \approx f'(a) + (x - a)f''(a)$$
$$f''(x) \approx f''(a) + (x - a)f'''(a)$$
Then,

$$
\begin{align*}
f(a + h) &\approx f(a) + \int_{a}^{a + h}{[f'(a) + (x - a)f''(a)]}\ dx \\
         &\approx f(a) + hf'(a) + \frac{h^2}{2}{f''(a)}
\end{align*}
$$

This process can be repeated to obtain higher order approximations to $f(a + h)$.

For the $(n - 1)$ th order approximation:

$$
f(a + h) \approx f(a) + hf'(a) + \frac{h^2}{2}{f''(a)}+...+\frac{h^{n-1}}{(n-1)!}f^{n-1}(a)
$$

The error associated with the approximation can be shown to be:
$$R_n(h) = \frac{h^n}{n!}f^{n}(\xi)$$
where $\xi \in [a,\ a+h]$

Then,
$$f(a + h) = f(a) + hf'(a) + \frac{h^2}{2}{f''(a)}+...+\frac{h^{n-1}}{(n-1)!}f^{n-1}(a) + R_n(h)$$

Substituting $x = a + h$,
$$f(x) = f(a) + (x-a)f'(a) + \frac{(x-a)^2}{2}{f''(a)}+...+\frac{(x-a)^{n-1}}{(n-1)!}f^{n-1}(a) + R_n(x-a)$$
