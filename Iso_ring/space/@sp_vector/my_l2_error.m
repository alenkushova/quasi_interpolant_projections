% MY_L2_ERROR: Evaluate the error in L^2 norm.
%
%   errl2 = my_l2_error (space, msh, u, uex, spex, coeff)
  
%
% INPUT:
%
%   space: object defining the space of discrete functions (see sp_scalar)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   u:     vector of dof weights
%   uex:   vector of dof weights of the exact solution
%   spex:  object defining the space of discrete functions (see sp_scalar)
%   coeff: handle function for the coeffitient of the exact solution 
%
% OUTPUT:
%
%     errl2:  error in L^2 norm
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011, 2015 Rafael Vazquez
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function errl2 = my_l2_error (space, msh, u, uex, spex, coeff)
   
  if (numel(u) ~= space.ndof)
    error ('Wrong size of the vector of degrees of freedom of u')
  end
  
  if (numel(uex) ~= spex.ndof)
    error ('Wrong size of the vector of degrees of freedom of uex')
  end
  
  errl2 = 0;
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'gradient', false);
    spex_col= sp_evaluate_col (spex, msh_col, 'value', true, 'gradient', false);
    
    if (nargin == 6)
      for idim = 1:msh.rdim
        x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
      end
      coeffs = coeff (x{:});
      coeffs = reshape(coeffs, 1, msh_col.nqn, msh_col.nel);
      coeffs = repmat(coeffs,msh_col.ndim,1,1);
    else
      coeffs = ones (msh_col.ndim, msh_col.nqn, msh_col.nel);
    end

    errl2 = errl2 + (my_l2_error (sp_col, msh_col, u, uex, spex_col, coeffs)).^2;
  end
  
  errl2 = sqrt (errl2);

end
