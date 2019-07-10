classdef CovarianceFunctionTest < matlab.unittest.TestCase
    %COVARIANCEFUNCTIONTEST Unit tests for the class
    %opencossan.common.inputs.covariancefunction
    % See also: OPENCOSSAN.COMMON.INPUTS.COVARIANCEFUNCTION
    
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
            Xcov = opencossan.common.inputs.stochasticprocess.CovarianceFunction();
            testCase.assertClass(Xcov,'opencossan.common.inputs.stochasticprocess.CovarianceFunction');
        end
        
        function constructorFullObject(testCase)
            Xcov = opencossan.common.inputs.stochasticprocess.CovarianceFunction('Format','structure',...
                'IsFunction',true,...
                'InputNames',{'t1','t2'},...
                'FullFileName',fullfile(opencossan.OpenCossan.getRoot(),'test','unit','data','common','inputs','CovarianceFunction','expcovfunction.m'),...
                'OutputNames',{'fcov'});
            testCase.assertEqual(Xcov.OutputNames,{'fcov'});
            testCase.assertEqual(Xcov.InputNames,{'t1' 't2'});
            testCase.assertTrue(Xcov.IsFunction);
        end
        
        function constructorShouldFailWithLessThanTwoInputs(testCase)
            testCase.assertError(@()opencossan.common.inputs.stochasticprocess.CovarianceFunction('Format','structure',...
                'IsFunction',true,...
                'InputNames',{'t1'},...
                'FullFileName',fullfile(opencossan.OpenCossan.getRoot(),'test','unit','data','common','inputs','CovarianceFunction','expcovfunction.m'),...
                'OutputNames',{'fcov'}),...
                'OpenCossan:CovarianceFunction:wrongNumberOfInputs');
        end
        
        function constructorShouldFailMoreLessThanTwoInputs(testCase)
            testCase.assertError(@()opencossan.common.inputs.stochasticprocess.CovarianceFunction('Format','structure',...
                'IsFunction',true,...
                'InputNames',{'t1' 't2' 't3'},...
                'FullFileName',fullfile(opencossan.OpenCossan.getRoot(),'test','unit','data','common','inputs','CovarianceFunction','expcovfunction.m'),...
                'OutputNames',{'fcov'}),...
                'OpenCossan:CovarianceFunction:wrongNumberOfInputs');
        end
        
        function constructorShouldFailWithMultipleOutputs(testCase)
            testCase.assertError(@()opencossan.common.inputs.stochasticprocess.CovarianceFunction('Format','structure',...
                'IsFunction',true,...
                'InputNames',{'t1', 't2'},...
                'FullFileName',fullfile(opencossan.OpenCossan.getRoot(),'test','unit','data','common','inputs','CovarianceFunction','expcovfunction.m'),...
                'OutputNames',{'fcov','more'}),...
                'OpenCossan:CovarianceFunction:wrongNumberOfOutputs');
        end
        
        %% evaluate
        function compute(testCase)
            Xcov = opencossan.common.inputs.stochasticprocess.CovarianceFunction('Format','structure',...
                'IsFunction',true,...
                'InputNames',{'t1','t2'},...
                'FullFileName',fullfile(opencossan.OpenCossan.getRoot(),'test','unit','data','common','inputs','CovarianceFunction','expcovfunction.m'),...
                'OutputNames',{'fcov'});
            
            MX = [0 0; 0 1; 1 0; 1 1];
            Vcov = Xcov.compute(MX);
            testCase.assertLength(Vcov,size(MX,2));
        end
        
        function computeMonoDimensionalSP(testCase)
            Xcov = opencossan.common.inputs.stochasticprocess.CovarianceFunction('Format','structure',...
                'IsFunction',true,...
                'InputNames',{'t1','t2'},...
                'FullFileName',fullfile(opencossan.OpenCossan.getRoot(),'test','unit','data','common','inputs','CovarianceFunction','expcovfunction.m'),...
                'OutputNames',{'fcov'});
            
            MX = [0 0; 1 1;];
            Vcov = Xcov.compute(MX);
            testCase.assertSize(Vcov,[2 1]);
        end
        
        function computeShouldFailForInvalidDimensions(testCase)
            Xcov = opencossan.common.inputs.stochasticprocess.CovarianceFunction('Format','structure',...
                'IsFunction',true,...
                'InputNames',{'t1','t2'},...
                'FullFileName',fullfile(opencossan.OpenCossan.getRoot(),'test','unit','data','common','inputs','CovarianceFunction','expcovfunction.m'),...
                'OutputNames',{'fcov'});
            MX = [0 0 1; 0 1 1; 1 0 1];
            testCase.assertError(@() Xcov.compute(MX),'OpenCossan:input:CovarianceFunction:evaluate');
        end
        
        function computeShouldFailWithoutInput(testCase)
            Xcov = opencossan.common.inputs.stochasticprocess.CovarianceFunction('Format','structure',...
                'IsFunction',true,...
                'InputNames',{'t1','t2'},...
                'FullFileName',fullfile(opencossan.OpenCossan.getRoot(),'test','unit','data','common','inputs','CovarianceFunction','expcovfunction.m'),...
                'OutputNames',{'fcov'});
            testCase.assertError(@() Xcov.compute(),'MATLAB:minrhs');
        end
    end
    
end
