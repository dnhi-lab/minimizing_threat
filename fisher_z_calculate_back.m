function [out] = fisher_z_calculate_back( in )

out = ( exp( 2*in ) - 1) ./ ( exp( 2*in ) + 1 );

end