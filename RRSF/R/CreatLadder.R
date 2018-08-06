CreatLadder <-
function( Ntotal, pRatio=0.75, Nmin=5 )
{
    x <- vector()
    x[1] <- Ntotal
    for( i in 1:100 )
    {
        pp <- round(x[i] * pRatio)
        if( pp == x[i] )
        {
            pp <- pp-1
        }          
        
        if( pp >= Nmin )
        {
            x[i+1] <- pp
        } else
        {
            break
        }
    }
    
    x
}
