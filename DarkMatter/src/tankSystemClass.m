classdef tankSystemClass < handle

    properties

        propellantTank
        pressurantTank

        m
        cg
        
        order
        type
    end

    methods

        function obj = tankSystemClass(inputTag, input)

            obj.addPropellantTank(inputTag, input);

            if input.(inputTag).isPressurized
                obj.addPressurantTank(input.(inputTag).pressurant, input);
            end
        end

        function addPropellantTank(obj, inputTag, input)
            if isempty(obj.propellantTank)
                obj.propellantTank = propellantTankClass(inputTag, input);
            else
                error('Propellant tank already defined!')
            end
        end

        function removePropellantTank(obj)
            if ~isempty(obj.propellantTank)
                obj.propellantTank = [];
            else
                error('No propellant tank is defined!')
            end
        end

        function addPressurantTank(obj, inputTag, input)
            if isempty(obj.pressurantTank)
                obj.pressurantTank = pressurantTankClass(inputTag, input);
            else
                error('Pressurant tank already defined!')
            end
        end

        function removePressurantTank(obj)
            if ~isempty(obj.pressurantTank)
                obj.pressurantTank = [];
            else
                error('No pressurant tank is defined!')
            end
        end

        function getCG(obj)





            if obj.isPressurized

                obj.pressurant.getCG();

                if strcmp(obj.pressurantOrder, 'fwd')

                    for i = 1:length(obj.mprop)

                        obj.cg(i) = (obj.pressurant.cg(i)*obj.pressurant.m(i) ...
                                  + (obj.pressurant.tank.l + obj.pressurant.offset + obj.cgprop(i))*obj.mprop(i)) ...
                                  / (obj.mprop(i) + obj.pressurant.m(i));
                    end

                elseif strcmp(obj.pressurantOrder, 'aft')

                    for i = 1:length(obj.mprop)

                        obj.cg(i) = ((obj.tank.l + obj.pressurant.offset + obj.pressurant.cg(i))*obj.pressurant.m(i) ...
                                  + (obj.cgprop(i)*obj.mprop(i))) ...
                                  / (obj.mprop(i) + obj.pressurant.m(i));
                    end
                end

                obj.m   = obj.mprop + obj.pressurant.m;

            else
                obj.cg  = obj.cgprop;
                obj.m   = obj.mprop;
            end
        end


    end

end
