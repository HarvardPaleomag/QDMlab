function chi = get_chi(final_fits)
    chi = final_fits.leftNeg.chiSquares+final_fits.leftPos.chiSquares+...
          final_fits.rightNeg.chiSquares+final_fits.rightPos.chiSquares;
end