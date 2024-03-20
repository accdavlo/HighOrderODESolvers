# HighOrderODESolvers
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

This is a repository for the course on "Topics in high order accurate time integration methods" held in the Doctoral School at SISSA Trieste ([Course Page 2024](https://www.math.sissa.it/course/phd-course/topics-high-order-accurate-time-integration-methods-0)). 


The code developed in the course is available in this repository under the MIT License, without any warranty.



Goal
======
The course aims to make the student aware of the cutting-edge algorithms used to numerically perform time integration and their properties in order to optimally choose the proper time integrator according to the considered phyisical model.

Prerequisites
======
* Basics of numerical analysis
* The course is based on Python Notebooks. They can be run locally on your laptop **or** online through Google Colab (**suggested**, you just need a Google account). Choose **one of the following possibilities**:
   * **Google Colab** installed on Google Drive (all online on Google servers) (**suggested** for unexperienced users) is an app that runs on Google Drive of your Google account. The installation is straightforward at the opening of the first Google Colab project. Try opening the [Google Colab link](https://colab.research.google.com/github/accdavlo/HighOrderODESolvers/blob/master/Chapter%200%20Test.ipynb) of Chapter 0 and run all the cells, follow the instructions that are inside the notebook. More info on Google Colab at [Google Colab Main Page](https://research.google.com/colaboratory) and [FAQ of Google Colab](https://research.google.com/colaboratory/faq.html).
   * **Jupyter Notebook on your laptop** (**unrecommended** if you do not use already python on your laptop and you are not a bit experienced with jupyter notebooks): you need to install python3 on your laptop, jupyter notebook and few modules (numpy, scipy, nodepy and matplotlib), best if installed with pip.
      * python3 : [Installation with Anaconda](https://www.anaconda.com/products/individual#Downloads)
      * jupyter notebook: [Installation with Anaconda and with pip](https://test-jupyter.readthedocs.io/en/latest/install.html)
      * Add the modules: numpy, scipy, nodepy and matplotlib: with pip ```pip install package-name``` with conda ```conda install package-name```.
      * Run the notebook [Guide](https://test-jupyter.readthedocs.io/en/latest/running.html#running)
      * Clone or download this repository (top right green button) 
      * Test the installation with the [Chapter 0 Notebook](Chapter%200%20Test.ipynb)


Schedule 2024
========
| Day                 | Time            |Room |
|  :------------      |:--------------- |:--- |
| Tuesday,   March 12 | 11:00 to 13:00  | 131 |
| Tuesday,   March 12 | 14:00 to 16:00  | 131 |
| Wednesday, March 13 | 11:00 to 13:00  | **132** |
| Wednesday, March 13 | 14:00 to 16:00  | 131 |
| Thursday,  March 14 | 14:00 to 18:00  | 131 |
| Tuesday,   March 19 | 11:00 to 13:00  | 131 |
| Wednesday, March 20 | 11:00 to 13:00  | **133** |
| Wednesday, March 20 | 14:00 to 16:00  | 131 |
| Thursday,  March 21 | 14:00 to 16:00  | 131 |


Microsoft Teams link
===========
The couse can be attended online at [this link](https://teams.microsoft.com/l/meetup-join/19%3ameeting_MmYxMDg5YjYtNTFmMy00ZmJhLWJmMDQtZmQyM2Y5NjA4MWRh%40thread.v2/0?context=%7b%22Tid%22%3a%22e4dd3336-ea1f-432c-b1e1-6966e8584f1b%22%2c%22Oid%22%3a%22cff80bd0-453f-4c4f-a6fd-61e7af40553d%22%7d).

Recordings
===========
* [Tuesday,    March 12, 11:00 to 13:00](https://sissa-my.sharepoint.com/:v:/g/personal/dtorlo_sissa_it/EaInWF2qW8VFsquajqZ8RHIBl8TryrNoNrTdPeFQ0E3aLw?e=zdjicg)
* [Tuesday,    March 12, 14:00 to 16:00](https://sissa-my.sharepoint.com/:v:/g/personal/dtorlo_sissa_it/EUTbdtcwpL9Pkn1NS4u4ZBEB4HM2ou1s4UKHX9b3PsFp8A?e=lQfGBt)
* [Wednesday,  March 13, 11:00 to 13:00](https://sissa-my.sharepoint.com/:v:/g/personal/dtorlo_sissa_it/ERBQe4wwPXtFizjBtzHmwcABom1rQZ_AxHn9LOF4EwD02g)
* [Wednesday,  March 13, 14:00 to 16:00](https://sissa-my.sharepoint.com/:v:/g/personal/dtorlo_sissa_it/EWvyiziA6gBBliTqViZENXsBPPYvFPkksuMEBW43TmqC0g)
* [Thursday,   March 14, 14:00 to 18:00](https://sissa-my.sharepoint.com/:v:/g/personal/dtorlo_sissa_it/EeceHltRGhhOjxODoQGtpFMB8uCLlFddnt4MKr2ONWmGng?e=spbs9K)
* [Tuesday,    March 19, 11:00 to 13:00](https://sissa-my.sharepoint.com/:v:/g/personal/dtorlo_sissa_it/Eb6f7KX2omlMnhqVn5taJ6cBAWcB4eqNbbaxCps_Ht3AkA?e=6wGbAp)
* [Wednesday,  March 20, 11:00 to 13:00](https://sissa-my.sharepoint.com/:v:/g/personal/dtorlo_sissa_it/Ea6L7I9TtbVKiooUXOFWip8BA_GnaSm9sJUSvbpGb0dlrA?e=aCYvDU)
* [Wednesday,  March 20, 14:00 to 16:00](https://sissa-my.sharepoint.com/:v:/g/personal/dtorlo_sissa_it/EZ7_Z_V8hrVFmEcunvRBvfIB76hEAkpNOQ5A69EM0Kyq9Q?e=bkRd0U)



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



## Disclaimer for the code

Everything is provided as-is and without warranty. Use at your own risk!
