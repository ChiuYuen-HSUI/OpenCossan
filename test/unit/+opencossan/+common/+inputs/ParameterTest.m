classdef ParameterTest < matlab.unittest.TestCase
    %PARAMETERTEST Unit tests for the class opencossan.common.inputs.Parameter
    % See also: OPENCOSSAN.COMMON.INPUTS.PARAMETER
    
    % This file is part of *OpenCossan*: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License.
    %
    % *OpenCossan* is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with openCOSSAN.  If not, see <http://www.gnu.org/licenses/>.
    % =====================================================================
    
    methods (Test)
        %% constructor
        function constructorEmpty(testCase)
            Xpar = opencossan.common.inputs.Parameter();
            testCase.verifyClass(Xpar, 'opencossan.common.inputs.Parameter');
            testCase.verifyNumElements(Xpar, 1);
            testCase.verifyLength(Xpar.Value, 0);
            testCase.verifyEqual(Xpar.Nelements, 0);
        end
        
        function constructorShouldSetDescription(testCase)
            Xpar = opencossan.common.inputs.Parameter('description','Description');
            testCase.verifyEqual(Xpar.Description,"Description");
        end
        
        function constructorShouldSetValue(testCase)
            Xpar = opencossan.common.inputs.Parameter('value', 5);
            testCase.verifyEqual(Xpar.Value, 5);
        end
        
        function constructorShouldValidateInput(testCase)
            % String validation for Description
            testCase.verifyError(@()opencossan.common.inputs.Parameter('description', cell(1)),...
                'MATLAB:UnableToConvert');
            testCase.verifyError(@()opencossan.common.inputs.Parameter('description', rand(2)),...
                'MATLAB:type:InvalidInputSize');
            % Numeric validation for Value
            testCase.verifyError(@()opencossan.common.inputs.Parameter('value', 'c'),...
                'MATLAB:validators:mustBeNumeric');
        end
        
        function constructorClassShouldNotInheritFromHandle(testCase)
            Xpar = opencossan.common.inputs.Parameter();
            testCase.verifyFalse(ishandle(Xpar));
        end
        
        
        %% display
        function checkDisplayWorks(testCase)
            % Check for single Object
            Xpar = opencossan.common.inputs.Parameter('description', 'Test Object',...
                'value', magic(4));
            testPhrase = [Xpar.Description;...
                num2str(Xpar.Nelements)];
            worksSingle = testOutput(Xpar,testPhrase);
            % Check for array of objects
            Xpar = [opencossan.common.inputs.Parameter(); ...
                opencossan.common.inputs.Parameter()];
            testPhrase = ["Parameter array with properties:";...
                "Description";...
                "Nelements";...
                "Value"];
            worksMulti = testOutput(Xpar,testPhrase);
            
            works = worksSingle && worksMulti;
            testCase.assertTrue(works);
        end
    end
end