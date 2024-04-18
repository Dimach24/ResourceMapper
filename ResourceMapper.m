classdef ResourceMapper<handle

    properties
        resourceGrid    % main grid
        mu              % Subcarrier spacing configuration
        frameCount      % amount of frames in resourceGrid
        isCycledPrefixExtended  % cycled prefix flag
    end

    methods
        
        function pbchDmRs=preparePbchDmRs(obj,pbchDmRs)
            % pre-mapping staff
            pbchDmRs=pbchDmRs+0; %todo
        end
        function pbch=preparePbch(obj,pbch)
             % pre-mapping staff
            pbch=pbch+0; %todo
        end
    
        function obj = ResourceMapper()
            % empty constructor;
        end

        
        function obj=createResourceGrid(obj,mu,frameCount,isCycledPrefixExtended)
            arguments
                obj
                mu                      (1,1)
                frameCount              (1,1)
                isCycledPrefixExtended  (1,1) =false
            end

            if(isCycledPrefixExtended) % ONLY FOR MU==2
                %   12 sym per slot, 10*4 slot per frame = 480 symb per frame
                obj.resourceGrid=zeros(4,480*frameCount); % FIXME 
            else
                obj.resourceGrid=zeros(240,2^mu*140*frameCount);% FIXME
            end
            % initializing
            obj.mu=mu;
            obj.isCycledPrefixExtended=isCycledPrefixExtended;
            obj.frameCount=frameCount;
        end
        
        %%% ===============================================================
      
        
        % SS MAPPING
        function addPssToResourceGrid(obj,PssSignal,offset)
        %adds PSS to resource matrix
            arguments
                obj 
                PssSignal (1,127)
                offset (1,1) = 2 %optional
            end
                obj.resourceGrid(57:183,1+offset) = fft(PssSignal.').';
        end
        
        function  addSssToResourceGrid(obj, SssSignal,offset)
        %adds SSS to resource matrix
            arguments
                obj
                SssSignal (1,127)
                offset (1,1) = 4 %optional
            end
            obj.resourceGrid(57:183,1+offset) = fft(SssSignal.').';
        end
        

        function addPbchToResourceGrid(obj,NCellId,pbch,offset)
            arguments
                obj
                NCellId
                pbch
                offset =0
            end

            % nu parameter for shift of DM-RS
            nu=mod(NCellId,4);
            % pre-mapping staff
            pbch=obj.preparePbch(pbch);
            
            %throwing out dmrs indexes
            indexes=find(mod(1:1:240,4)~=(nu+1));

            % mapping first 180 PBCH
            obj.resourceGrid(indexes,1+offset)=pbch(1:180);
            % mapping last 180
            obj.resourceGrid(indexes,3+offset)=pbch(end-179:end);
            % mapping arround SSS
            indexes=indexes(indexes<49 | indexes>192);
            obj.resourceGrid(indexes,2+offset)=pbch(181:181+71);
        end


        %PBCH DM-RS MAPPING
        function addPbchDmRsToResourceGrid(obj,NCellId,pbchDmRs,offset)
            arguments
                obj
                NCellId
                pbchDmRs
                offset =0
            end
            % nu parameter for shift of DM-RS
            nu=mod(NCellId,4);
            % pre-mapping staff
            pbchDmRs=obj.preparePbchDmRs(pbchDmRs);
            % first 60 PBCH DM-RS
            dmrs=pbchDmRs(1:60);
            
            % indexes array
            indexes=find(mod(1:1:240,4)==1);
            % mapping 1st part
            obj.resourceGrid(indexes+nu,1+offset)=dmrs;
            % d---d---d---d … d---d---
            
            % last dm-rs block
            dmrs=pbchDmRs(end-59:end);
            % mapping 2nd part
            obj.resourceGrid(indexes+nu,3+offset)=dmrs;
            % d---d---d---d … d---d---

            % next dm-rs block (24 elements are splitted into two blocks)
            dmrs=pbchDmRs(62:62+23);
            indexes=indexes(indexes<46 | indexes>192); % throwing SSS area
            % mapping 3rd part
            obj.resourceGrid(indexes+nu,2+offset)=dmrs;
            % d---d---d--…-SSS-…-d---d--d
        end
    end
    
end