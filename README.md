# slsqp4j

Slsqp4j is a Java wrapper around the popular `SLSQP` nonlinear optimizer included in SciPy. The API mimics SciPy's in order to ease the 
translation of problems from Python to the JVM. 

The bulk of the solving is done in `slsqp.f90` which was written by Dieter Kraft and described in <a href="#ref1">[1]</a> 
& <a href="#ref2">[2]</a>.

## Building
Building Slsqp4j depends on both gcc and gfortran. 

### Ubuntu
You can install both with the command `sudo apt install gcc gfortran`. Additionally, your `JAVA_HOME`  must point to your JDK install directory. 

To build Slsqp4j, simply run `gradle clean build` in the project root directory.

### Mac OSX
This procedure has been tested with adoptopenjdk8 ( `brew tap AdoptOpenJDK/openjdk && brew install --cask adoptopenjdk8` ) on BigSur.

You can install both with the command `brew install gradle gcc gfortran`. 

Additionally, your `JAVA_HOME`  must point to your JDK real install directory (ie a path like
`/Library/Java/JavaVirtualMachines/adoptopenjdk-8.jdk/Contents/Home/` and not a path like
`/Library/Internet Plug-Ins/JavaAppletPlugin.plugin/Contents/Home` that doesn't contain the `include` library. 

To build Slsqp4j, simply run `gradle clean build` in the project root directory.

**Errors & solutions:**
`/usr/local/Cellar/gcc/10.2.0/lib/gcc/10/gcc/x86_64-apple-darwin19/10.2.0/include-fixed/stdio.h:78:10: fatal error: _stdio.h: No such file or directory
   78 | #include <_stdio.h>
      |          ^~~~~~~~~~
compilation terminated.`
==> run `xcode-select --install`
    Then run XCode and install additional components

### Windows

(_Note: The following steps were performed using [scoop command-line installer](https://www.scoop.sh)._)

You can install both with the command `scoop install mingw-winlibs`. Additionally, your `JAVA_HOME`  must point to your JDK install directory.

To build Slsqp4j, simply run `gradle clean build` in the project root directory.

## Usage

Create an objective function that implements the `Vector2ScalarFunc` interface:
```Java
    public static class ObjectiveFunction implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            // for example
            return x[0] * x[1];
        }
    }
```

Specify one or more constraints:
```Java
    public static final class VectorConstraintFunction implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {x[0] - x[1]};
        }
    }

    final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
        .withConstraintType(ConstraintType.EQ)
        .withConstraintFunction(new VectorConstraintFunction())
        .build();
```
 
To perform the optimization, you must construct an instance of an `Slsqp` object. You do this using the Builder pattern:
```Java
final Slsqp slsqp = new Slsqp.SlsqpBuilder()
    .withObjectiveFunction(new ObjectiveFunction())
    .addVectorConstraint(constraint)
    .build();
```
Then simply call `optimize` passing in an initial guess vector:
```Java
final OptimizeResult result = slsqp.minimize(new double[] {1, -1});
```

The returned `OptimizeResult` contains information about the state of the solver. If `result.success()` returns true,
the solver is complete and the vector contained in `result.resultVec()` is the point at which the function is minimized.

Below is a comparison showing the complete usage of Slsqp4j's API vs. SciPy's `optimize` API.
<table>
<th>
Slsqp4j
</th>
</table>

```Java
final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
    .withConstraintType(ConstraintType.EQ)
    .withConstraintFunction((x, arg) -> x[0] - x[1])
    .build();
final Slsqp slsqp = new Slsqp.SlsqpBuilder()
    .withObjectiveFunction((x, arg) -> x[0] * x[1])
    .addVectorConstraint(constraint)
    .build();
final OptimizeResult result = slsqp.minimize(new double[]{1, -1});
```

<table>
<th>
SciPy
</th>
</table>

```python
res = minimize(lambda d: d[0] * d[1], [1, -1], method='SLSQP', 
        constraints={'type': 'ineq', 'fun': lambda x: x[0] - x[1]})
```


The API is slightly more verbose than SciPy's one due to Java's type safety, but the similarities should hopefully be apparent. 
For more usage examples refer to the tests in [SlsqpTests.java](./slsqp4j/src/test/java/com/skew/slsqp4j/SlsqpTests.java). For a complete list of the `SlsqpBuilder` 
parameters consult the documentation in [Slsqp.java](./slsqp4j/src/main/java/com/skew/slsqp4j/Slsqp.java).

### Thread Safety
Since the `SLSQP` algorithm is iterative, it is assumed that an instance of `Slsqp` will not be shared among threads, thus instances are *not* thread-safe. Rather, an instance of `Slsqp` should be constructed once, with the parameters of the optimization problem given to the builder, and then repeated calls to `slsqp.optimize()` should be made until a value of `true` is returned on a call to `success()` on the returned `OptimizeResult` instance.

## Tests
The majority of the tests in [SlsqpTests.java](./slsqp4j/src/test/java/com/skew/slsqp4j/SlsqpTests.java) were ported from SciPy's SLSQP tests.

## License
Slsqp4j is released under the [BSD license](https://github.com/skew-markets/slsqp4j/blob/master/LICENSE.txt).

## References
<ol>
<li id="ref1">Dieter Kraft, "A software package for sequential quadratic
programming", Technical Report DFVLR-FB 88-28, Institut für
Dynamik der Flugsysteme, Oberpfaffenhofen, July 1988.</li>

<li id="ref2">Dieter Kraft, "Algorithm 733: TOMP–Fortran modules for optimal
control calculations," ACM Transactions on Mathematical Software,
vol. 20, no. 3, pp. 262-281 (1994).</li>
</ol>
