# ImproveNLS
Efforts to improve the functioning of the R `nls()` function for nonlinear least 
squares estimation. This is part of the ***Google Summer of Code*** program for 2021.

### Aim of the Project
Our aim is to focus attention on actions that can be realized in future versions of the installable R system.
Given the very wide range of applications and features of the `nls()` function, however, it is almost certain that only some possibilities will be addressed in the project. 

### Project Title
**Improvements to nls()**

### Project Mentors
- John C. Nash, Telfer School of Management, University of Ottawa
- Heather Turner, University of Warwick
- Hans W. Borchers

### Project Student 
- Arkajyoti Bhattacharjee

### Project Outcomes
- **Final Reports:**

  1. [RefactoringNLS](https://gitlab.com/nashjc/improvenls/-/blob/master/RefactoringNLS.pdf)
  <br/>

- **Other (informal) Reports:**<br/>

  2. [PkgFromRbase](https://gitlab.com/nashjc/improvenls/-/blob/master/PkgFromRbase.pdf): explanation of the construction of the [nlspkg](https://gitlab.com/nashjc/improvenls/-/tree/master/nlspkg) from the code in R-base.
  3. [DerivsNLS](https://gitlab.com/nashjc/improvenls/-/blob/master/DerivsNLS.pdf): document to explain different ways in which Jacobian information is supplied to nonlinear least squares computation in R.
  4. [MachineSummary](https://gitlab.com/nashjc/improvenls/-/blob/master/MachineSummary.pdf): informal investigation of ways to report the characteristics and identity of machines running tests.
  5. [VarietyInNonlinearLeastSquaresCodes](https://gitlab.com/nashjc/improvenls/-/blob/master/VarietyInNonlinearLeastSquaresCodes.pdf): review of the different algorithms and the many choices in their implementation for nonlinear least squares.
  6. [ImproveNLS.bib](https://gitlab.com/nashjc/improvenls/-/blob/master/ImproveNLS.bib): consolidated BibTex bibliography for all documents in this project, possibly with wider application to nonlinear least squares in general
  7. [WorkingDocument4ImproveNLS](https://gitlab.com/nashjc/improvenls/-/blob/master/WorkingDocument4ImproveNLS.pdf): project diary
  <br/>

- **Codes (and documentation):**<br/>

  8. [nlsj](https://gitlab.com/nashjc/improvenls/-/tree/master/nlsj): A refactoring of the `nls()` functionality
  9. [nlsCompare](): an R package to compare existing and new packages' functions for nonlinear least squares 
  10. [nlspkg](https://gitlab.com/nashjc/improvenls/-/tree/master/nlspkg): a packaged version of the `nls()` code from R-base.
  11. [nlsalt](https://gitlab.com/nashjc/improvenls/-/tree/master/nlsalt): attempt to mirror `nls()` behaviour **entirely in R**
  12. [nls-changes-for-small-residuals-in-nls-R-4.0.2.zip](https://gitlab.com/nashjc/improvenls/-/blob/master/nls-changes-for-small-residuals-in-nls-R-4.0.2.zip): collected material for the fix by JN to the **relative
offset convergence criterion failure** when there are small residuals in problems sent to `nls()`. 
  13. [nlsralt](https://github.com/ArkaB-DS/nlsCompare): a modified version of Nash and Murdoch package `nlsr` with improvements discovered as a result of this project.

### Contribute
We welcome further collaboration on this project. 

