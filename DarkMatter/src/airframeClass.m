classdef    airframeClass < handle

    properties (Access = public)
        input
        design
        utilities
        aero
    end

    methods
        function obj    = airframeClass(input)

            % Populate general properties
            obj.input           = input;                            % Input
            obj.aero.drag_file  = 'Drag_Data_Houbolt_Jr.csv';       %drag file to be read
            obj.drag_model();
        end

        function drag_model(obj)
            drag_data       = csvread(obj.aero.drag_file);
            obj.aero.cd     = fit(drag_data(:,1),drag_data(:,2),'cubicinterp');
        end
    end
end