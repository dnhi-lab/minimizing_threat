function [out] = fisher_z_calculate( in )

out = 0.5 * log( (1 + in ) ./ ( 1 - in ) );

end