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
        function obj = addPssToResourceGrid(obj,PssSignal,offset)
        %adds PSS to resource matrix
            arguments
                obj 
                PssSignal (1,127)
                offset (1,1) = 2 %optional
            end
            for i = 0:obj.frameCount*10*2^obj.mu-1 %frame amount * 10 subframes * 2^mu slots
                % one slot has 2 PSS
                obj.resourceGrid(57:183,1+offset+14*i) = fft(PssSignal.').';
                obj.resourceGrid(57:183,1+offset+6+14*i) = fft(PssSignal.').';
            end
        end
        
        function  addSssToResourceGrid(obj, SssSignal,offset)
        %adds SSS to resource matrix
            arguments
                obj
                SssSignal (1,127)
                offset (1,1) = 4 %optional
            end
            %frame amount * 10 subframes * 2^mu slots
            for i = 0:obj.frameCount*10*2^obj.mu-1 
                % one slot has 2 SSS
                obj.resourceGrid(57:183,1+offset+14*i) = fft(SssSignal.').';
                obj.resourceGrid(57:183,1+offset+6+14*i) = fft(SssSignal.').';
            end
        end
        
        % PBCH MAPPING
        function obj = addPbchToResourceGrid(obj,NCellId,pbch)
           % for each slot in frames
           for i=0:1:obj.frameCount*2^obj.mu*10-1
               % 2 pbch parts per slot
                obj=obj.addPbchPartToResourceGrid(NCellId,pbch(i+1,:),2+i*14);
                obj=obj.addPbchPartToResourceGrid(NCellId,pbch(i+1,:),8+i*14);
           end 
        end

        function obj = addPbchPartToResourceGrid(obj,NCellId,pbch,shift)
            arguments
                obj
                NCellId
                pbch
                shift =0
            end

            % nu parameter for shift of DM-RS
            nu=mod(NCellId,4);
            % pre-mapping staff
            pbch=obj.preparePbch(pbch);
            
            %throwing out dmrs indexes
            indexes=find(mod(1:1:240,4)~=(nu+1));

            % mapping first 180 PBCH
            obj.resourceGrid(indexes,2+shift)=pbch(1:180);
            % mapping last 180
            obj.resourceGrid(indexes,4+shift)=pbch(end-179:end);
            % mapping arround SSS
            indexes=indexes(indexes<49 | indexes>192);
            obj.resourceGrid(indexes,3+shift)=pbch(181:181+71);
        end

        function obj=addPbchDmRsToResourceGrid(obj,NCellId,pbchDmRs)
            % for each slot in grid
            for i=0:1:obj.frameCount*2^obj.mu*10-1
                 % twice in each slot of the subframe
                obj=obj.addPbchDmRsPartToResourceGrid(NCellId,pbchDmRs(i+1,:),2+i*14);
                obj=obj.addPbchDmRsPartToResourceGrid(NCellId,pbchDmRs(i+1,:),8+i*14);
            end
        end

        %PBCH DM-RS MAPPING
        function obj=addPbchDmRsPartToResourceGrid(obj,NCellId,pbchDmRs,shift)
            arguments
                obj
                NCellId
                pbchDmRs
                shift =0
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
            obj.resourceGrid(indexes+nu,2+shift)=dmrs;
            % d---d---d---d … d---d---
            
            % last dm-rs block
            dmrs=pbchDmRs(end-59:end);
            % mapping 2nd part
            obj.resourceGrid(indexes+nu,4+shift)=dmrs;
            % d---d---d---d … d---d---

            % next dm-rs block (24 elements are splitted into two blocks)
            dmrs=pbchDmRs(62:62+23);
            indexes=indexes(indexes<46 | indexes>192); % throwing SSS area
            % mapping 3rd part
            obj.resourceGrid(indexes+nu,3+shift)=dmrs;
            % d---d---d--…-SSS-…-d---d--d
        end
    end
    
end