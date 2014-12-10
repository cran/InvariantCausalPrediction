pvalfunc <-
function( x, y) 2*min( t.test(x,y)$p.value, var.test(x,y)$p.value)
