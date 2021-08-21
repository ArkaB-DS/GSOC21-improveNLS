# ImproveNLS
Efforts to improve the functioning of the R `nls()` function for nonlinear least 
squares estimation. This is part of the ***`Google Summer of Code`*** program for 2021.

## Aim of the Project
Our aim is to focus attention on actions that can be realized in future versions of the installable R system.
Given the very wide range of applications and features of the `nls()` function, however, it is almost certain that only some possibilities will be addressed in the project. 

## Project Title
**Improvements to nls()**

## Project Organization
**The R Project for Statistical Computing**

## Project Mentors
- John C. Nash, Telfer School of Management, University of Ottawa
- Heather Turner, University of Warwick

## Project Student 
- Arkajyoti Bhattacharjee

## Project Outcomes
- **Final Reports:**
 
  1. [RefactoringNLS](https://gitlab.com/nashjc/improvenls/-/blob/master/RefactoringNLS.pdf)
  2. [nlsCompareArticle](https://github.com/ArkaB-DS/GSOC21-improveNLS/blob/master/nlsCompareArticle/nlsCompareArticle.pdf)

- **Other (informal) Reports:**
 
  3. [PkgFromRbase](https://github.com/ArkaB-DS/GSOC21-improveNLS/blob/master/PkgFromRbase.pdf): explanation of the construction of the [nlspkg](https://gitlab.com/nashjc/improvenls/-/tree/master/nlspkg) from the code in R-base.
  4. [DerivsNLS](https://github.com/ArkaB-DS/GSOC21-improveNLS/blob/master/DerivsNLS.pdf): document to explain different ways in which Jacobian information is supplied to nonlinear least squares computation in R.
  5. [MachineSummary](https://github.com/ArkaB-DS/GSOC21-improveNLS/blob/master/MachineSummary.pdf): informal investigation of ways to report the characteristics and identity of machines running tests.
  6. [VarietyInNonlinearLeastSquaresCodes](https://github.com/ArkaB-DS/GSOC21-improveNLS/blob/master/VarietyInNonlinearLeastSquaresCodes.pdf): review of the different algorithms and the many choices in their implementation for nonlinear least squares.
  7. [ImproveNLS.bib](https://github.com/ArkaB-DS/GSOC21-improveNLS/blob/master/ImproveNLS.bib): consolidated BibTex bibliography for all documents in this project, possibly with wider application to nonlinear least squares in general
  8. [WorkingDocument4ImproveNLS](https://github.com/ArkaB-DS/GSOC21-improveNLS/blob/master/WorkingDocument4ImproveNLS.pdf): project diary

- **Codes (and documentation):**
  
  9. [nlsj](https://github.com/ArkaB-DS/nlsj): A refactoring of the `nls()` functionality
  10. [nlsCompare](https://github.com/ArkaB-DS/nlsCompare): an R package to compare existing and new packages' functions for nonlinear least squares 
  11. [nlspkg](https://github.com/ArkaB-DS/GSOC21-improveNLS/tree/master/nlspkg): a packaged version of the `nls()` code from R-base.
  12. [nlsalt](https://github.com/ArkaB-DS/GSOC21-improveNLS/tree/master/nlsalt): attempt to mirror `nls()` behaviour **entirely in R**
  13. [nls-changes-for-small-residuals-in-nls-R-4.0.2.zip](https://github.com/ArkaB-DS/GSOC21-improveNLS/blob/master/nls-changes-for-small-residuals-in-nls-R-4.0.2.zip): collected material for the fix by JN to the **relative
offset convergence criterion failure** when there are small residuals in problems sent to `nls()`. 
  14. [nlsralt](https://github.com/ArkaB-DS/GSOC21-improveNLS/tree/master/nlsralt): a modified version of Nash and Murdoch package `nlsr` with improvements discovered as a result of this project.

## Mentions

We would like to thank -
- **Hans Werner Borchers** for the contributions to the proposal of the project

