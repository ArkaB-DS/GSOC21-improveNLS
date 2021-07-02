## "FIXME": move to 'base' and make .Internal or even .Primitive
setNames <- function(object = nm, nm)
{
    names(object) <- nm
    object
}
