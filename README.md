# slsqp4j

Slsqp4j is a Java wrapper around the popular `SLSQP` optimizer included in Scipy. The bulk of the solving
is done by the code written in `slsqp.f90` which was originally written by Dieter Kraft and described in <a href="#ref1">[1]</a> 
& <a href="#ref2">[2]</a>.

### Building
To build slsqp4j you will need to install both gcc and gfortran. 
This can be done with the command 

`sudo apt install gcc gfortran`

Additionally, your `JAVA_HOME`  must point to your JDK install directory. 

### Usage
To use Slslqp4j you must construct an instance of an `Slsqp4j` object. You do this using the Builder pattern:
```
final Slsqp slsqp = new Slsqp.SlsqpBuilder()
    .withObjectiveFunction(new Fun())
    .withJacobian(new Jac())
    .addScalarConstraint(constraint)
    .build();
```

The builder accepts an objective function, as well as a Jacobian, some constraints, bounds, error tolerance, and a maximum number
of iterations for the Slsqp solver to perform. Some of these parameters are optional. The Slsqp solver can 
solve unconstrained and unbounded problems.


Below is a side-by-side comparison showing Slsqp4j's API vs. Scipy's `optimize` api.
<table>
<tr>
<th>
Slsqp4j
</th>
<th>
Scipy
</th>
</tr>

<tr>
<td>
<pre>
final VectorConstraint constraint1 = new VectorConstraint.VectorConstraintBuilder()
    .withConstraintType(ConstraintType.EQ)
    .withConstraintFunction(new Fecon2())
    .build();
final VectorConstraint constraint2 = new VectorConstraint.VectorConstraintBuilder()
    .withConstraintType(ConstraintType.EQ)
    .withConstraintFunction(new Fecon())
    .build();
final VectorConstraint constraint3 = new VectorConstraint.VectorConstraintBuilder()
    .withConstraintType(ConstraintType.INEQ)
    .withConstraintFunction(new Fieqcon2())
    .build();<br>
final double[] lowerBounds = new double[] {-0.8, -1};
final double[] upperBounds = new double[] {1, 0.8};
final double[][] bounds = new double[][] {lowerBounds, upperBounds};<br>
final Slsqp slsqp = new Slsqp.SlsqpBuilder()
    .withObjectiveFunction(new Fun(), -1)
    .withJacobian(new Jac())
    .withBounds(bounds)
    .addVectorConstraint(constraint1)
    .addVectorConstraint(constraint2)
    .addVectorConstraint(constraint3)
    .withAccuracy(defaultTol)
    .withMaxIterations(defaultMaxIter)
    .build();
final OptimizeResult result = slsqp.minimize(new double[]{-1.4, 0.9});
</pre>
</td>
<td>
<pre>
constraints = [
    {'type': 'eq', 'fun': self.f_eqcon2},
    {'type': 'eq', 'fun': self.f_eqcon},
    {'type': 'ineq', 'fun': self.f_ieqcon2},
] 
res = minimize(self.fun, [-1.4, 0.9], method='SLSQP',
       jac=self.jac, args=(-1.0, ),
       bounds=[(-0.8, 1.), (-1, 0.8)])

</pre>
</td>
</tr>
</table>

The API is slightly more verbose than the Scipy one due to Java's type safety, however the similarities should be apparent. 

### License
Slsqp4j is released under the [BSD license](https://github.com/skew-markets/slsqp4j/blob/master/LICENSE.txt).

### References
<ol>
<li id="ref1">Dieter Kraft, "A software package for sequential quadratic
programming", Technical Report DFVLR-FB 88-28, Institut für
Dynamik der Flugsysteme, Oberpfaffenhofen, July 1988.</li>

<li id="ref2">Dieter Kraft, "Algorithm 733: TOMP–Fortran modules for optimal
control calculations," ACM Transactions on Mathematical Software,
vol. 20, no. 3, pp. 262-281 (1994).</li>
</ol>
