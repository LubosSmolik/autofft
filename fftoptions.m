classdef fftoptions
    properties
        % Property validation uses an enumeration class. This technique
        % (i) provides inexact, case-insensitive matching for unambiguous
        % char vectors or string scalars, and (ii) converts inexact matches
        % to correct values

        FFTLength (1,1)
        Mode   (1,1) propMode
        Window (1,1) propWindow
        Dome (1,1) string {mustBeMember(Dome, ["red", "blue"])} = "red"
    end

    methods
        % Constructor
        function obj = fftoptions(args)
            arguments
                % Ensures that 
                args.Mode   string
                args.Window {string, double}

                % Allow the constructor to accept an argument args that is
                % an instance of the fftoptions class
                args.?fftoptions
            end

            for prop = string(fieldnames(args))'
                obj.(prop) = args.(prop);
            end
        end

        % Ensure that strings are returned rather than enumerations
        function out = get.Mode(obj)
            out = string(obj.Mode);
        end

        function out = get.Window(obj)
            out = string(obj.Window);
        end
    end
end
