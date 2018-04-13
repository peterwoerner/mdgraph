function [str1 ] = datesave( str )
%Creates a string with the date behind it for saving files
%
c = clock;
str1 = strcat(str,'_', num2str(c(1)), num2str(c(2)), num2str(c(3)), num2str(c(4)), num2str(c(5)))

end

