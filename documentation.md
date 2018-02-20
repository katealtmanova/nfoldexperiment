---


---

<h1 id="nfoldip-technical-documentation">NFoldIP technical documentation</h1>
<p>We provide a class to handle instances of n-fold integer programming.<br>
Most configuration is done via the arguments of the <code>__init__</code> function (see below).<br>
The majority of computational work is done in the <code>_find_good_step</code> function, which searches for a good augmenting step of a specified length.<br>
The algorithm computes such a step for each step-length from a set of step-lengths Gamma and then uses the best of these to augment the current solution; it terminates when no further augmenting steps can be found.<br>
The <code>_find_good_step</code> procedure is implemented using the dynamic programming approach of [1]; because of its poor performance it can also be solved as an ILP subproblem, which is solved using a MILP solver such as GLPK, Coin-OR or Gurobi; see the function <code>solve</code>.</p>
<p>All of the important implemented functions are described below.</p>
<p>The code was written in Python/Cython as a SageMath’s module using the package 4ti2 — a software package for algebraic, geometric and combinatorial problems on linear spaces. Available at <a href="http://www.4ti2.de">www.4ti2.de</a>.</p>
<h2 id="implementation---functions">IMPLEMENTATION - functions:</h2>
<p>The <code>__init__</code> function - takes these arguments:</p>
<ol>
<li><code>self</code></li>
<li><code>A</code> - diagonal matrix (vectors of machines when scheduling)</li>
<li><code>D</code> - upper matrix (matrix of jobs when scheduling)</li>
<li><code>n</code> - (number of jobs when scheduling)</li>
<li><code>b</code> -right hand-side vector</li>
<li><code>l</code> - lower bound</li>
<li><code>u</code> - upper bound</li>
<li><code>w</code> - objective function</li>
<li><code>verbose</code> - the identifier of logging, standard system of<br>
loggers - <code>notset</code>, <code>debug</code>, <code>info</code>, <code>warning</code>, <code>error</code> (default), <code>critical</code> (see<br>
<a href="https://docs.python.org/2/library/logging.html">https://docs.python.org/2/library/logging.html</a> ), optional.</li>
<li><code>graver_complexity</code> - possible values are: <code>“exact”</code> - computes the value of Graver complexity from the definition, <code>“approximate”</code> (default) - computes an upper bound using a formula, integer - user-input; turns the algorithm into a heuristic</li>
<li><code>current_solution</code> - optional, if no current solution is given, the algorithm computes the initial solution using an auxiliary instance,</li>
<li><code>experimental</code> - optional, <code>False</code> (default) uses dynamic programming (slow), <code>True</code> uses a MILP solver to solve a subinstance equivalent to the DP, <code>"ng1"</code> uses a MILP solver to find augmenting steps with 1-norm bounded by <code>graver_complexity</code>, <code>"nginfty"</code> uses a MILP solver to find augmenting steps with \infty-norm bounded by <code>graver_complexity</code>.</li>
<li><code>instancename</code> - optional, a string which will be used for output file logging, defaults to <code>"instancename"</code></li>
<li><code>gamma</code> - optional, a string, <code>"best"</code> uses the “best step” augmentation strategy, <code>"logarithmic"</code> (default) uses the “approximate best step” augmentation strategy, <code>"unit"</code> uses the “any step” augmentation strategy,</li>
<li><code>solver</code> - optional, a string, determines which MILP solver to use if <code>experimental is not False</code>, <code>"GLPK"</code> (default), <code>"Coin"</code>, <code>"Gurobi"</code>.</li>
</ol>
<blockquote>
<p><strong>Implementation:</strong></p>
</blockquote>
<blockquote>
<ul>
<li>checking the validity of the given data by a function check_validity_of_data, following must hold:</li>
</ul>
<blockquote>
<p>general size of matrices:<br>
<img src="https://lh3.googleusercontent.com/qvX3mVPqOhcma_LnSqK4uZ-4olQnBYPfmI9X3TP7JBMmFyCNAHzjqJRqt7i7LzPGLi1u8XvFFXZd" alt="enter image description here" title="general size of matrices"><br>
together with the sizes of l/u bounds, x and b:<img src="https://lh3.googleusercontent.com/dAUwdfJjT0pZwB1_Hk_XngEUHih6wT4i1BBg5oWTiX5pJvAOXeFQMfFApVp-OaGmCUhViqbvmMvc" alt="enter image description here" title="picture of sizes"></p>
</blockquote>
</blockquote>
<blockquote>
<ul>
<li>lower and upper bounds can’t be infinity</li>
<li>and for every number l in the lower bound vector and for every number u in the upper bound vector must hold that l&lt;=u (for l,u in the same position)</li>
<li>the logging level + format, self arguments such as self.A, self.D, self.t, self.r, self.s etc. are being initialized</li>
<li>if any initial feasible solution was given, it is being set to self.current_solution; if there isn’t any, the initial feasible solution will be computed by the find_init_feasible_solution function later</li>
<li>graver complexity is  being computed and set to self.graver_complexity by a function approximate_graver_complexity, or exact_graver_complexity, or only the value is being applied (according to the graver_complexity option)</li>
<li>ZE to self.ZE is being set through the function <code>_construct_ZE</code> (only if <code>experimental == false</code>)</li>
<li>checking the validity of the data (whether the sizes of matrices, l/u bounds etc. make sense</li>
</ul>
</blockquote>
<p>Two functions for computing the graver complexity of given matrix (at most one of these is chosen according to the tenth (graver_complexity) argument in the <strong>init</strong> function):</p>
<p><strong>exact_graver_complexity</strong> returns an exact value of the Graver complexity of the given data from the definition</p>
<blockquote>
<p><strong>Implementation</strong>:</p>
</blockquote>
<blockquote>
<ul>
<li>takes the maximum entry of a matrix which  is product of multiplication of D by the graver basis of the matrix A ( the graver base has been computed using 4ti2/set as an integer in the init function)</li>
</ul>
</blockquote>
<p><strong>approximate_graver_complexity</strong> return an approximate value of the Graver complexity of the given matrices</p>
<blockquote>
<p><strong>Implementation</strong>:</p>
</blockquote>
<blockquote>
<ul>
<li>computes graver complexity according to the following formula [2]<br>
<span class="katex--display"><span class="katex-display"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>g</mi><mo>(</mo><mi>A</mi><mo>)</mo><mo>≤</mo><mi>p</mi><mo>(</mo><mi>r</mi><mo>∗</mo><mi mathvariant="normal">∣</mi><mi mathvariant="normal">∣</mi><mi>D</mi><mo>∗</mo><mi>G</mi><mi>A</mi><mi mathvariant="normal">∣</mi><msub><mi mathvariant="normal">∣</mi><mi mathvariant="normal">∞</mi></msub><msup><mo>)</mo><mi>r</mi></msup></mrow><annotation encoding="application/x-tex">g(A)\leq p(r*|| D*GA||_{\infty})^{r} </annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.75em;"></span><span class="strut bottom" style="height: 1em; vertical-align: -0.25em;"></span><span class="base"><span style="margin-right: 0.03588em;" class="mord mathit">g</span><span class="mopen">(</span><span class="mord mathit">A</span><span class="mclose">)</span><span class="mrel">≤</span><span class="mord mathit">p</span><span class="mopen">(</span><span style="margin-right: 0.02778em;" class="mord mathit">r</span><span class="mbin">∗</span><span class="mord mathrm">∣</span><span class="mord mathrm">∣</span><span style="margin-right: 0.02778em;" class="mord mathit">D</span><span class="mbin">∗</span><span class="mord mathit">G</span><span class="mord mathit">A</span><span class="mord mathrm">∣</span><span class="mord"><span class="mord mathrm">∣</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: 0em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathrm mtight">∞</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mclose"><span class="mclose">)</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height: 0.714392em;"><span class="" style="top: -3.113em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span style="margin-right: 0.02778em;" class="mord mathit mtight">r</span></span></span></span></span></span></span></span></span></span></span></span></span></span><br>
where:</li>
</ul>
<blockquote>
<ul>
<li>g(A) is the graver  complexity of matrix A</li>
<li>p is the number of elements in the graver basis of A</li>
<li>GA is the graver basis of A</li>
</ul>
</blockquote>
</blockquote>
<p>Function for computing Z(E) - <strong>construct_ZE</strong>. Z(E) is the sum of at most Graver complexity elements of the matrix A. This function is also called from the <code>__init__</code> function.</p>
<blockquote>
<p><strong>Implementation:</strong></p>
<ul>
<li>at first creates vector of zeros (it is definitely in ZE)</li>
<li>then in a three inner for cycles happen following:</li>
</ul>
<blockquote>
<ul>
<li>1st cycle:   graver complexity times new empty set is created</li>
<li>2nd cycle:   depending on the size of yet computed unique elements of ZE</li>
<li>3rd cycle:   two vectors are computed (graver complexity times) — it’s a sum/difference of one vector from yet computed ZE with a vector from graver basis of <code>A</code></li>
</ul>
</blockquote>
</blockquote>
<p>Finding feasible solution - <code>find_init_feasible_solution</code>. It computes the initial solution if it has not been given in the <code>__init__</code> function. It consists of two methods - <strong>create_auxiliary_program</strong>, which creates the instance of an auxiliary program. Then there is a method for solving the aux instance - <strong>solve_auxiliary_program</strong>.</p>
<blockquote>
<p><strong>Implementation  (</strong><code>create_auxiliary_program</code><strong>)</strong>:</p>
</blockquote>
<blockquote>
<ul>
<li>at first it constructs two matrices:</li>
</ul>
<blockquote>
<p><span class="katex--display"><span class="katex-display"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msub><mi>D</mi><mn>2</mn></msub><mo>=</mo><mo>[</mo><mi>D</mi><msub><mi>Z</mi><mrow><mi>s</mi><mi>r</mi></mrow></msub><msub><mi>Z</mi><mrow><mi>s</mi><mi>r</mi></mrow></msub><msub><mi>I</mi><mi>s</mi></msub><mo>−</mo><msub><mi>I</mi><mi>s</mi></msub><mo>]</mo></mrow><annotation encoding="application/x-tex">D_{2} = [D Z_{sr} Z_{sr} I_{s} -I_{s}]</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.75em;"></span><span class="strut bottom" style="height: 1em; vertical-align: -0.25em;"></span><span class="base"><span class="mord"><span style="margin-right: 0.02778em;" class="mord mathit">D</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.301108em;"><span class="" style="top: -2.55em; margin-left: -0.02778em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathrm mtight">2</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mrel">=</span><span class="mopen">[</span><span style="margin-right: 0.02778em;" class="mord mathit">D</span><span class="mord"><span style="margin-right: 0.07153em;" class="mord mathit">Z</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: -0.07153em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathit mtight">s</span><span style="margin-right: 0.02778em;" class="mord mathit mtight">r</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mord"><span style="margin-right: 0.07153em;" class="mord mathit">Z</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: -0.07153em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathit mtight">s</span><span style="margin-right: 0.02778em;" class="mord mathit mtight">r</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mord"><span style="margin-right: 0.07847em;" class="mord mathit">I</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: -0.07847em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathit mtight">s</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mbin">−</span><span class="mord"><span style="margin-right: 0.07847em;" class="mord mathit">I</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: -0.07847em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathit mtight">s</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mclose">]</span></span></span></span></span></span><br>
<span class="katex--display"><span class="katex-display"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msub><mi>A</mi><mn>2</mn></msub><mo>=</mo><mo>[</mo><mi>A</mi><msub><mi>I</mi><mi>r</mi></msub><mo>−</mo><msub><mi>I</mi><mi>r</mi></msub><msub><mi>Z</mi><mrow><mi>r</mi><mi>s</mi></mrow></msub><msub><mi>Z</mi><mrow><mi>r</mi><mi>s</mi></mrow></msub><mo>]</mo></mrow><annotation encoding="application/x-tex">A_{2} = [A I_{r} -I_{r} Z_{rs} Z_{rs}]</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.75em;"></span><span class="strut bottom" style="height: 1em; vertical-align: -0.25em;"></span><span class="base"><span class="mord"><span class="mord mathit">A</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.301108em;"><span class="" style="top: -2.55em; margin-left: 0em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathrm mtight">2</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mrel">=</span><span class="mopen">[</span><span class="mord mathit">A</span><span class="mord"><span style="margin-right: 0.07847em;" class="mord mathit">I</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: -0.07847em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span style="margin-right: 0.02778em;" class="mord mathit mtight">r</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mbin">−</span><span class="mord"><span style="margin-right: 0.07847em;" class="mord mathit">I</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: -0.07847em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span style="margin-right: 0.02778em;" class="mord mathit mtight">r</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mord"><span style="margin-right: 0.07153em;" class="mord mathit">Z</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: -0.07153em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span style="margin-right: 0.02778em;" class="mord mathit mtight">r</span><span class="mord mathit mtight">s</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mord"><span style="margin-right: 0.07153em;" class="mord mathit">Z</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.151392em;"><span class="" style="top: -2.55em; margin-left: -0.07153em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span style="margin-right: 0.02778em;" class="mord mathit mtight">r</span><span class="mord mathit mtight">s</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span><span class="mclose">]</span></span></span></span></span></span></p>
<blockquote>
<ul>
<li>where <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msub><mi>M</mi><mrow><mi>a</mi><mi>b</mi></mrow></msub></mrow><annotation encoding="application/x-tex">M_{ab}</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.68333em;"></span><span class="strut bottom" style="height: 0.83333em; vertical-align: -0.15em;"></span><span class="base"><span class="mord"><span style="margin-right: 0.10903em;" class="mord mathit">M</span><span class="msupsub"><span class="vlist-t vlist-t2"><span class="vlist-r"><span class="vlist" style="height: 0.336108em;"><span class="" style="top: -2.55em; margin-left: -0.10903em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span class="mord mathit mtight">a</span><span class="mord mathit mtight">b</span></span></span></span></span><span class="vlist-s">​</span></span><span class="vlist-r"><span class="vlist" style="height: 0.15em;"></span></span></span></span></span></span></span></span></span> means a matrix of <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>a</mi></mrow><annotation encoding="application/x-tex">a</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.43056em;"></span><span class="strut bottom" style="height: 0.43056em; vertical-align: 0em;"></span><span class="base"><span class="mord mathit">a</span></span></span></span></span> rows and <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>b</mi></mrow><annotation encoding="application/x-tex">b</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.69444em;"></span><span class="strut bottom" style="height: 0.69444em; vertical-align: 0em;"></span><span class="base"><span class="mord mathit">b</span></span></span></span></span> columns, <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>I</mi></mrow><annotation encoding="application/x-tex">I</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.68333em;"></span><span class="strut bottom" style="height: 0.68333em; vertical-align: 0em;"></span><span class="base"><span style="margin-right: 0.07847em;" class="mord mathit">I</span></span></span></span></span> is an identity matrix, <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>Z</mi></mrow><annotation encoding="application/x-tex">Z</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.68333em;"></span><span class="strut bottom" style="height: 0.68333em; vertical-align: 0em;"></span><span class="base"><span style="margin-right: 0.07153em;" class="mord mathit">Z</span></span></span></span></span> is a matrix full of zeros</li>
</ul>
</blockquote>
</blockquote>
<ul>
<li>then then it computes new lower and upper bounds:</li>
</ul>
<blockquote>
<ul>
<li>lower bound vector consists of <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mo>(</mo><mi>t</mi><mo>+</mo><mn>2</mn><mo>∗</mo><mi>r</mi><mo>+</mo><mn>2</mn><mo>∗</mo><mi>s</mi><mo>)</mo><mo>∗</mo><mi>n</mi></mrow><annotation encoding="application/x-tex">(t+2*r+2*s)*n</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.75em;"></span><span class="strut bottom" style="height: 1em; vertical-align: -0.25em;"></span><span class="base"><span class="mopen">(</span><span class="mord mathit">t</span><span class="mbin">+</span><span class="mord mathrm">2</span><span class="mbin">∗</span><span style="margin-right: 0.02778em;" class="mord mathit">r</span><span class="mbin">+</span><span class="mord mathrm">2</span><span class="mbin">∗</span><span class="mord mathit">s</span><span class="mclose">)</span><span class="mbin">∗</span><span class="mord mathit">n</span></span></span></span></span> zeros</li>
<li>upper bound vector consist of <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mo>(</mo><mi>t</mi><mo>+</mo><mn>2</mn><mo>∗</mo><mi>r</mi><mo>+</mo><mn>2</mn><mo>∗</mo><mi>s</mi><mo>)</mo><mo>∗</mo><mi>n</mi></mrow><annotation encoding="application/x-tex">(t+2*r+2*s)*n</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.75em;"></span><span class="strut bottom" style="height: 1em; vertical-align: -0.25em;"></span><span class="base"><span class="mopen">(</span><span class="mord mathit">t</span><span class="mbin">+</span><span class="mord mathrm">2</span><span class="mbin">∗</span><span style="margin-right: 0.02778em;" class="mord mathit">r</span><span class="mbin">+</span><span class="mord mathrm">2</span><span class="mbin">∗</span><span class="mord mathit">s</span><span class="mclose">)</span><span class="mbin">∗</span><span class="mord mathit">n</span></span></span></span></span> times the max value in the <code>self.b</code> vector</li>
</ul>
</blockquote>
<ul>
<li>creates new objective function:</li>
</ul>
<blockquote>
<ul>
<li>it consist of a vector of <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>t</mi></mrow><annotation encoding="application/x-tex">t</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.61508em;"></span><span class="strut bottom" style="height: 0.61508em; vertical-align: 0em;"></span><span class="base"><span class="mord mathit">t</span></span></span></span></span> times zero and <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mn>2</mn><mo>∗</mo><mi>r</mi><mo>+</mo><mn>2</mn><mo>∗</mo><mi>s</mi></mrow><annotation encoding="application/x-tex">2*r+2*s</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.64444em;"></span><span class="strut bottom" style="height: 0.72777em; vertical-align: -0.08333em;"></span><span class="base"><span class="mord mathrm">2</span><span class="mbin">∗</span><span style="margin-right: 0.02778em;" class="mord mathit">r</span><span class="mbin">+</span><span class="mord mathrm">2</span><span class="mbin">∗</span><span class="mord mathit">s</span></span></span></span></span> times one which is n-times copied</li>
</ul>
</blockquote>
<ul>
<li>makes the initial feasible solution of the auxiliary program</li>
</ul>
<blockquote>
<ul>
<li>the vector of the initial feasible solution consists of 2 types of vectors:</li>
</ul>
<blockquote>
<ul>
<li>the first type has <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>t</mi><mo>+</mo><mn>2</mn><mi>r</mi><mo>+</mo><mn>2</mn><mi>s</mi></mrow><annotation encoding="application/x-tex">t+2r+2s</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.64444em;"></span><span class="strut bottom" style="height: 0.72777em; vertical-align: -0.08333em;"></span><span class="base"><span class="mord mathit">t</span><span class="mbin">+</span><span class="mord mathrm">2</span><span style="margin-right: 0.02778em;" class="mord mathit">r</span><span class="mbin">+</span><span class="mord mathrm">2</span><span class="mord mathit">s</span></span></span></span></span> numbers and each number is a value from the lower vector (on the corresponding position) if it’s not minus infinity, elif it’s a value from the upper bound vector if it’s not an infinity, elif it is zero; then is this vector filled with positive numbers form the <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>b</mi><mo>[</mo><mn>0</mn><mo>]</mo></mrow><annotation encoding="application/x-tex">b[0]</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.75em;"></span><span class="strut bottom" style="height: 1em; vertical-align: -0.25em;"></span><span class="base"><span class="mord mathit">b</span><span class="mopen">[</span><span class="mord mathrm">0</span><span class="mclose">]</span></span></span></span></span> vector (or zeros when not positive), then the negative numbers from the <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>b</mi><mo>[</mo><mn>0</mn><mo>]</mo></mrow><annotation encoding="application/x-tex">b[0]</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.75em;"></span><span class="strut bottom" style="height: 1em; vertical-align: -0.25em;"></span><span class="base"><span class="mord mathit">b</span><span class="mopen">[</span><span class="mord mathrm">0</span><span class="mclose">]</span></span></span></span></span> vector (or zeros when not negative) and the same with the <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>b</mi><mo>[</mo><mn>1</mn><mo>]</mo></mrow><annotation encoding="application/x-tex">b[1]</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.75em;"></span><span class="strut bottom" style="height: 1em; vertical-align: -0.25em;"></span><span class="base"><span class="mord mathit">b</span><span class="mopen">[</span><span class="mord mathrm">1</span><span class="mclose">]</span></span></span></span></span> vector (positive and negative part)</li>
<li>the second type of the vector of the init feasible solution is actually the same as the first type but the first part is different, there are only zeros<br>
-note: the vector <code>b</code> consists of two parts - one corresponds to the <code>D</code> matrix, the second part corresponds to the <code>A</code> matrix</li>
</ul>
</blockquote>
</blockquote>
<ul>
<li>finally it returns an instance of NFoldIP</li>
</ul>
</blockquote>
<blockquote>
<p><strong>Implementation (</strong><code>solve_auxiliary_program</code><strong>):</strong></p>
<ul>
<li>its argument is an auxiliary instance of NFoldIP with its initial solution</li>
<li>uses the algorithm to minimize the auxiliary variables in order to have the initial solution for the main program</li>
<li>then it checks whether the auxiliary vars are zero — if yes, we have an initial feasible solution (returns the init feasible solution), otherwise the main program has no feasible solution (returns None)</li>
</ul>
</blockquote>
<p>If the initial feasible solution exists we are searching for the augmenting steps by the function <strong>find_graverbest_step</strong>.</p>
<blockquote>
<p><strong>Implementation:</strong></p>
<ul>
<li>at the beginning the generator of Gamma of gammas has to be computed (see the function <code>construct_Gamma below</code>)</li>
<li>then for the gamma in Gamma:</li>
</ul>
<blockquote>
<p>while the dot product of vector <code>w</code> with the <code>good_step</code> (which was computed from the <code>_find_good_step(gamma)</code> function) is not greater or equal to zero:</p>
<blockquote>
<ul>
<li>try to take bigger gamma in order to prolong the good step</li>
</ul>
</blockquote>
</blockquote>
<ul>
<li>at the end the function returns the maximum step, which is <code>gamma*good_step</code> with the best dot product with corresponding <code>w</code></li>
</ul>
</blockquote>
<p>Now there is a short look into the construction of Gamma - <strong>construct_Gamma</strong>. It is an iterator of gammas for possible extension of the lengths of the feasible steps. It is used in the <code>find_graverbest_step</code> function as it has been mentioned before.</p>
<blockquote>
<p><strong>Implementation:</strong></p>
</blockquote>
<blockquote>
<ul>
<li>it is an iterator of logarithmic values</li>
<li>it makes a vector of the difference from upper-lower bounds, chooses the biggest element of the difference, the max value which this iterator returns is the number <code>n</code>(floor) for which holds following:</li>
</ul>
<blockquote>
<p><span class="katex--display"><span class="katex-display"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><msup><mi>n</mi><mi>r</mi></msup><mo>=</mo><mi>m</mi><mi>a</mi><mi>x</mi><mi>v</mi><mi>a</mi><mi>l</mi><mi>u</mi><mi>e</mi></mrow><annotation encoding="application/x-tex">n^{r} = max value</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.714392em;"></span><span class="strut bottom" style="height: 0.714392em; vertical-align: 0em;"></span><span class="base"><span class="mord"><span class="mord mathit">n</span><span class="msupsub"><span class="vlist-t"><span class="vlist-r"><span class="vlist" style="height: 0.714392em;"><span class="" style="top: -3.113em; margin-right: 0.05em;"><span class="pstrut" style="height: 2.7em;"></span><span class="sizing reset-size6 size3 mtight"><span class="mord mtight"><span style="margin-right: 0.02778em;" class="mord mathit mtight">r</span></span></span></span></span></span></span></span></span><span class="mrel">=</span><span class="mord mathit">m</span><span class="mord mathit">a</span><span class="mord mathit">x</span><span style="margin-right: 0.03588em;" class="mord mathit">v</span><span class="mord mathit">a</span><span style="margin-right: 0.01968em;" class="mord mathit">l</span><span class="mord mathit">u</span><span class="mord mathit">e</span></span></span></span></span></span></p>
</blockquote>
<ul>
<li>the first value is 1, every next value is two times the previous value if it is lower than the <code>maxvalue</code>, otherwise <code>StopIteration</code> is raised</li>
</ul>
</blockquote>
<p>There is the option not to use <code>find_good_step(gamma)</code>, but its experimental version - <strong>find_good_step_ecperimental(gamma, time limit)</strong>.</p>
<blockquote>
<p><strong>Implementation:</strong></p>
<ul>
<li>TODO</li>
</ul>
</blockquote>
<p>And finally the crucial function for solving a given n-fold program: <code>solve(“solver”)</code>. It takes one argument which says which solver should be used.</p>
<ul>
<li>if <code>"native"</code>, it uses the <code>find_graverbest_step</code> to compute the solution</li>
<li>if <code>“GLPK”</code> it creates a MILP instance and solves it with GLPK</li>
</ul>
<p>If we choose the native solver the function <strong>native_solve</strong> is called. It solves the given problem with the implemented algorithm.</p>
<blockquote>
<p><strong>Implementation:</strong></p>
<ul>
<li>it finds and applies the graver best steps by the function <code>find_graverbest_step</code> until there are no more graver best steps</li>
</ul>
</blockquote>
<p>The second option for solving calls the function <strong>glpk_solve</strong>. This function builds MILP in a standard form and uses GLPK for solving the problem.</p>
<blockquote>
<p><strong>Implementation:</strong></p>
<ul>
<li>sets <code>maximization = False</code>, add new nonnegative variable <code>d</code></li>
<li>builds the form of <span class="katex--inline"><span class="katex"><span class="katex-mathml"><math><semantics><mrow><mi>A</mi><mi>d</mi><mo>=</mo><mi>b</mi><mo separator="true">,</mo><mi>d</mi><mo>≥</mo><mn>0</mn></mrow><annotation encoding="application/x-tex">Ad=b, d≥0</annotation></semantics></math></span><span class="katex-html" aria-hidden="true"><span class="strut" style="height: 0.69444em;"></span><span class="strut bottom" style="height: 0.88888em; vertical-align: -0.19444em;"></span><span class="base"><span class="mord mathit">A</span><span class="mord mathit">d</span><span class="mrel">=</span><span class="mord mathit">b</span><span class="mpunct">,</span><span class="mord mathit">d</span><span class="mrel">≥</span><span class="mord mathrm">0</span></span></span></span></span></li>
<li>and minimizes <code>wd</code></li>
<li>it does not respect any initial solution if given at the beginning</li>
</ul>
</blockquote>
<h2 id="sources">SOURCES:</h2>
<p>[1] R. Hemmecke, S. Onn and L. Romanchuk, “N-fold integer programming in cubic time,” Mathematical Programming, pp. 1–17, 2013.</p>
<p>[2] S. Onn, Nonlinear discrete optimization, Zurich Lectures in Advanced Mathematics, European Mathematical Society, 2010.</p>
<!--stackedit_data:&#10;eyJoaXN0b3J5IjpbNDIwODU3NjkxXX0=&#10;-->

