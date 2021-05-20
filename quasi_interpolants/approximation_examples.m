function [] = approximation_examples ()

iopt = 1;
while (iopt > 0 )
    fprintf (1, ...
           ['2D Lee-Lyche-Morken approximation examples menu:\n', ...
            '----------------------\n', ...
            '\n', ...
            '   (1) Example of approximation with Dirichlet boundary. \n \n', ...
            '   (2) Example of approximation with periodic boundary. \n \n']);
  
    iopt = input ('Please choose a number from above or press <Enter> to return: ');
    clc;
    if (iopt == 1)
        example_Dirichlet_Approx
    elseif (iopt == 2)
        example_Periodic_Approx
    end
end