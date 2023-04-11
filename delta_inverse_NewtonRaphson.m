function SolutionToNonlinearAlgebraicEquations = delta_inverse
% File  delta_inverse.m

global   qA qB qR;
global   ABSERR;

ReadUserInput;
VAR(1) = qA;
VAR(2) = qB;
VAR(3) = qR;
SolutionToNonlinearAlgebraicEquations = SolveSetOfNonlinearAlgebraicEquations(VAR',ABSERR);
PrintUserOutput;
ErrorInConvergenceOfSolution = CalculateFunctionEvaluatedAtX( SolutionToNonlinearAlgebraicEquations )
CloseOutputFilesAndTerminate;



%===========================================================================
function ReadUserInput
global   D Iny la lo lr x y z;
global   qA qB qR;
global   DEGtoRAD RADtoDEG COEF RHS ABSERR;

%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
D                               =  20;                    % UNITS               Constant
Iny                             =  10;                    % UNITS               Constant
la                              =  30;                    % UNITS               Constant
lo                              =  20;                    % UNITS               Constant
lr                              =  30;                    % UNITS               Constant
x                               =  27.5979;                    % UNITS               Constant
y                               =  43.4405;                    % UNITS               Constant
z                               =  33.1500;                    % UNITS               Constant

qA                              =  0.0;                    % UNITS               Guess
qB                              =  0.0;                    % UNITS               Guess
qR                              =  0.0;                    % UNITS               Guess

ABSERR                          =  1.0E-08;                %                     Absolute Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
Pi       = 3.141592653589793;
DEGtoRAD = Pi/180.0;
RADtoDEG = 180.0/Pi;

% Reserve space and initialize matrices
COEF = zeros(3,3);
RHS = zeros(1,3);

% Evaluate constants



%===========================================================================
function SolutionToNonlinearAlgebraicEquations = SolveSetOfNonlinearAlgebraicEquations( X, AbsoluteError )
alphaBest = 0.5;
numberOfGradientCalculations = 0;
while 1
   errorInSetOfFunctions = CalculateErrorInSetOfFunctions( X );
   if( errorInSetOfFunctions < AbsoluteError ) break; end
   numberOfGradientCalculations = numberOfGradientCalculations + 1;
   dX = CalculateGradientAtX( X );
   if( max(abs(dX)) == 0.0 ) warning('Nonlinear solver failed'); break;  end
   alpha = sort( [0;  0.8*alphaBest; 1.2*alphaBest; 1] );
   alphaBest = CalculateBestValueOfAlpha( X, alpha, dX );
   if( alphaBest == 0.0 ) warning('Nonlinear solver failed'); break;  end
   X = X + alphaBest*dX;
end
SolutionToNonlinearAlgebraicEquations = X;
numberOfGradientCalculations = numberOfGradientCalculations;



%===========================================================================
function alphaBest = CalculateBestValueOfAlpha( X, alpha, dX )
error(1) = CalculateErrorInSetOfFunctions( X + alpha(1)*dX );
error(2) = CalculateErrorInSetOfFunctions( X + alpha(2)*dX );
error(3) = CalculateErrorInSetOfFunctions( X + alpha(3)*dX );
error(4) = CalculateErrorInSetOfFunctions( X + alpha(4)*dX );
numberOfFunctionCalculations = 4;
maxDistanceBeforeReturning = 1.0/(2+length(X))^2;
while 1
   if( error(1) <= error(2) )
      if( numberOfFunctionCalculations > 20 ) alphaBest=0.0;  return;  end
      alpha(4)=alpha(3);  alpha(3)=alpha(2);  alpha(2)=0.5*(alpha(1)+alpha(2));
      error(4)=error(3);  error(3)=error(2);  error(2)=CalculateErrorInSetOfFunctions( X+alpha(2)*dX );
   elseif( error(2) <= error(3) )
      alpha(4)=alpha(3);  alpha(3)=alpha(2);  alpha(2)=0.5*(alpha(1)+alpha(2));
      error(4)=error(3);  error(3)=error(2);  error(2)=CalculateErrorInSetOfFunctions( X+alpha(2)*dX );
   elseif( error(3) <= error(4) )
      if( alpha(4)-alpha(1) <= maxDistanceBeforeReturning ) alphaBest=alpha(3); break; end
      space23=alpha(3)-alpha(2);  space34=alpha(4)-alpha(3);
      if( space23 >= space34 )
         alpha(1)=alpha(2);  alpha(2)=alpha(2)+0.5*space23;
         error(1)=error(2);  error(2)=CalculateErrorInSetOfFunctions( X+alpha(2)*dX );
      else
         alpha(1)=alpha(2);  alpha(2)=alpha(3);  alpha(3)=alpha(3)+0.5*space34;
         error(1)=error(2);  error(2)=error(3);  error(3)=CalculateErrorInSetOfFunctions( X+alpha(3)*dX );
      end
   else
      alphaBest=alpha(4);   break;
   end
   numberOfFunctionCalculations = numberOfFunctionCalculations + 1;
end
numberOfFunctionCalculations = numberOfFunctionCalculations;



%===========================================================================
function ErrorInSetOfFunctions = CalculateErrorInSetOfFunctions( X )
ErrorInSetOfFunctions = norm( CalculateFunctionEvaluatedAtX(X) );



%===========================================================================
function RHS = CalculateFunctionEvaluatedAtX( VAR )
global   D Iny la lo lr x y z;
global   qA qB qR;
global   RHS;

% Update variables with new values
qA = VAR(1);
qB = VAR(2);
qR = VAR(3);
RHS(1) = z + lr*(1-cos(qB)) - Iny - la*sin(qA) - lr*sin(qR);
RHS(2) = 0.5*D + y - 0.5*lo - la*cos(qA) - lr*cos(qR);
RHS(3) = x - lr*sin(qB);



%===========================================================================
function SolutionToAxEqualsB = CalculateGradientAtX( VAR )
global   la  lr;
global   qA qB qR
global   COEF RHS SolutionToAxEqualsB ;

CalculateFunctionEvaluatedAtX( VAR );
COEF(1,1) = la*cos(qA);
COEF(1,2) = -lr*sin(qB);
COEF(1,3) = lr*cos(qR);
COEF(2,1) = -la*sin(qA);
COEF(2,2) = 0;
COEF(2,3) = -lr*sin(qR);
COEF(3,1) = 0;
COEF(3,2) = lr*cos(qB);
COEF(3,3) = 0;

SolutionToAxEqualsB = pinv(COEF)*RHS';

%===========================================================================
function Output = PrintUserOutput

%===========================================================================
function WriteOutput( fileIdentifier, Output )
numberOfOutputQuantities = length( Output );
if numberOfOutputQuantities > 0
  for i=1:numberOfOutputQuantities
    fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
  end
  fprintf( fileIdentifier, '\n' );
end

%===========================================================================
function CloseOutputFilesAndTerminate
fprintf( 1, '\n The quantities that are returned to the calling function are:\n' );
fprintf( 1, ' qA (UNITS)    qB (UNITS)    qR (UNITS)   \n\n' );
