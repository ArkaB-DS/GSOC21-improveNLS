#     int i;
#     #define CONV_INFO_MSG(_STR_, _I_)				\
#     ConvInfoMsg(_STR_, i, _I_, fac, minFac, maxIter, convNew)
#  ?? how does this ConvInfoMsg work
#     
#     ?? Seems this is defining the behaviour i.e., a function definition
#     ?? so the code is NOT executed here, but when NON_CONV_FINIS?
#     ??  are called. 
#     #define NON_CONV_FINIS(_ID_, _MSG_)		\
#     if(warnOnly) {				\
#         warning(_MSG_);				\
#         return CONV_INFO_MSG(_MSG_, _ID_);      \
#     }						\
#     else					\
#     error(_MSG_);
# ?? here are my replacements
##      cat("ctrl$warnOnly =",ctrl$warnOnly,"\n")
##      if (ctrl$warnOnly){
##         warning("-- need ConvInfoMsg here ?? --")
##         return(NULL) #?? need to return proper info to match nls()
##      } else {
##         stop("--need a suitable msg to stop --??")
##      }
#     
#     #define NON_CONV_FINIS_1(_ID_, _MSG_, _A1_)	\
#     if(warnOnly) {				\
#         char msgbuf[1000];			\
#         warning(_MSG_, _A1_);			\
#         snprintf(msgbuf, 1000, _MSG_, _A1_);	\
#         return CONV_INFO_MSG(msgbuf, _ID_);	\
#     }						\
#     else					\
#     error(_MSG_, _A1_);
#     
#     #define NON_CONV_FINIS_2(_ID_, _MSG_, _A1_, _A2_)	\
#     if(warnOnly) {					\
#         char msgbuf[1000];				\
#         warning(_MSG_, _A1_, _A2_);			\
#         snprintf(msgbuf, 1000, _MSG_, _A1_, _A2_);	\
#         return CONV_INFO_MSG(msgbuf, _ID_);		\
#     }							\
#     else						\
#     error(_MSG_, _A1_, _A2_);
#======= End of defining the messages about non-convergence ======     
