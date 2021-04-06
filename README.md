# HighOrderODESolvers
This is a repository for the course on "High order accurate time integration methods"

[Course Page](https://www.adum.fr/script/formations.pl?mod=292601&site=edmi)

Goal
======
The course aims to make the student aware of the cutting-edge algorithms used to numerically perform time integration and their properties in order to optimally choose the proper time integrator according to the considered phyisical model.

Prerequisites
======
* Basics of numerical analysis
* Python Jupyter Notebook on your laptop (suggested if you already use python) or Google Colab installed on Google Drive (all online) (better for unexperienced users)
   * Jupyter Notebook: you need to install python3 on your laptop, jupyter notebook and few modules (numpy, scipy, nodepy and matplotlib), best if installed with pip.
      * python3 : [Installation with Anaconda](https://www.anaconda.com/products/individual#Downloads)
      * jupyter notebook: [Installation with Anaconda and with pip](https://test-jupyter.readthedocs.io/en/latest/install.html)
      * Add the modules: numpy, scipy, nodepy and matplotlib (if with pip you can run the first cell of the installation) with conda, type ```conda install package-name```.
      * Run the notebook [Guide](https://test-jupyter.readthedocs.io/en/latest/running.html#running)
      * Clone or download this repository (top right green button) 
      * Test the installation with the [Chapter 0 Notebook](Chapter%200%20Test.ipynb)
   * Google Colab is app that runs on Google Drive. After you installed it you can open the link of Google Colab below on your drive! [Google Colab Main Page](https://research.google.com/colaboratory), [FAQ of Google Colab](https://research.google.com/colaboratory/faq.html). Open the [Google Colab link](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%200%20Test.ipynb) of Chapter 0 and run all the cells.


Schedule
========
 * Session 1: 14/04/2021 9:00-12:00
 * Session 2: 14/04/2021 13:30-16:30
 * Session 3: 15/04/2021 9:00-12:00
 * Session 4: 15/04/2021 13:30-16:30

Program
======
 * Chapter 0: Test of the installation [Notebook](Chapter%200%20Test.ipynb) or [Google Colab](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%200%20Test.ipynb)
 * Chapter 1: Theory of ODEs: examples, existence and uniqueness of solutions. [Notebook](Chapter&#32;1&#32;Theory&#32;of&#32;ODE.ipynb) or [Google Colab](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%201%20Theory%20of%20ODE.ipynb)
 * Chapter 2: Explicit and implicit Euler, convergence, stability analysis, properties. [Notebook](Chapter&#32;2&#32;Classical&#32;Euler&#32;Methods.ipynb) or [Google Colab](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%202%20Classical%20Euler%20Methods.ipynb)
 * Chapter 3: High order classical methods [Notebook](Chapter&#32;3&#32;Classical&#32;High&#32;Order&#32;Methods.ipynb) or [Google Colab](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%203%20Classical%20High%20Order%20Methods.ipynb)
   * Runge--Kutta methods: construction, explicit, implicit, error and stability analysis, properties, the Butcher tableau.
   * Multistep methods: explicit, implicit methods, error and stability analysis, convergence.
 * Chapter 4: Entropy/energy conservative/dissipative high order schemes: relaxation Runge--Kutta methods. [Notebook](Chapter&#32;4&#32;Relaxation&#32;Runge--Kutta.ipynb) or [Google Colab](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%204%20Relaxation%20Runge--Kutta.ipynb)
 * Chapter 5: Iterative explicit arbitrarily high order methods. Deferred Correction (DeC), Arbitrary Derivative (ADER) methods, properties, stability, convergence analysis, similarities. [Notebook](Chapter&#32;5&#32;DeC&#32;and&#32;ADER.ipynb) or [Google Colab](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%205%20DeC%20and%20ADER.ipynb) and [Slides on ADER and DeC](Chapter5/latexSlides/ADERDeC_chapter5.pdf) 
 * Chapter 6: Unconditionally positivity preserving schemes and strong stability preserving schemes. Modified Patankar methods and strong strability preserving Runge Kutta and multistep methods.
[Notebook](Chapter&#32;6&#32;Positivity&#32;preserving&#32;schemes.ipynb) or [Google Colab](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%206%20Positivity%20preserving%20schemes.ipynb) and [Slides for modified Patankar](Chapter6/latexSlides/mPDeC_Chapter6.pdf)

Solutions of exercises
------
They are in folder [solutions](/solutions)
You can open them in Google Colab through
1. [Google Colab solutions](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/solutions/Chapter%201%20Theory%20of%20ODE.ipynb)
1. [Google Colab solutions](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/solutions/Chapter%202%20Classical%20Euler%20Methods.ipynb)
1. [Google Colab solutions](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/solutions/Chapter%203%20Classical%20High%20Order%20Methods.ipynb)
1. [Google Colab solutions](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/solutions/Chapter%204%20Relaxation%20Runge--Kutta.ipynb)
1. [Google Colab solutions](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/solutions/Chapter%205%20DeC%20and%20ADER.ipynb)
1. [Google Colab solutions](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/solutions/Chapter%206%20Positivity%20preserving%20schemes.ipynb)
