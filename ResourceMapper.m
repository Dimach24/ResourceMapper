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
      

        function addSsBlockToResourceGrid(obj,nCellId,pssSignal,sssSignal,pbch,pbchDmRs,t_offset,f_offset)
            arguments
                obj ResourceMapper
                nCellId int32
                pssSignal
                sssSignal
                pbch
                pbchDmRs
                t_offset=0
                f_offset=0
            end
            obj.addPssToResourceGrid(pssSignal,t_offset,f_offset);
            obj.addSssToResourceGrid(sssSignal,t_offset+2,f_offset);
            obj.addPbchDmRsToResourceGrid(nCellId,pbchDmRs,t_offset+1,f_offset);
            obj.addPbchToResourceGrid(nCellId,pbch,t_offset+1,f_offset)
        end

        % SS MAPPING
        function addPssToResourceGrid(obj,PssSignal,t_offset,f_offset)
        %adds PSS to resource matrix
            arguments
                obj 
                PssSignal (1,127)
                t_offset (1,1) = 0 %optional
                f_offset (1,1) = 0 %optional
            end
                obj.resourceGrid((57:183)+f_offset,1+t_offset) = fft(PssSignal.').';
        end
        
        function  addSssToResourceGrid(obj, SssSignal,t_offset,f_offset)
        %adds SSS to resource matrix
            arguments
                obj
                SssSignal (1,127)
                t_offset (1,1) = 0 %optional
                f_offset (1,1) = 0 %optional
            end
            obj.resourceGrid((57:183)+f_offset,1+t_offset) = fft(SssSignal.').';
        end
        
        function addPbchToResourceGrid(obj,NCellId,pbch,t_offset,f_offset)
            arguments
                obj
                NCellId
                pbch
                t_offset = 0
                f_offset = 0
            end

            % nu parameter for shift of DM-RS
            nu=mod(NCellId,4);
            % pre-mapping staff
            pbch=obj.preparePbch(pbch);
            
            %throwing out dmrs indexes
            indexes=find(mod(1:1:240,4)~=(nu+1));

            % mapping first 180 PBCH
            obj.resourceGrid(indexes+f_offset,1+t_offset)=pbch(1:180);
            % mapping last 180
            obj.resourceGrid(indexes+f_offset,3+t_offset)=pbch(end-179:end);
            % mapping arround SSS
            indexes=indexes(indexes<49 | indexes>192);
            obj.resourceGrid(indexes+f_offset,2+t_offset)=pbch(181:181+71);
        end

        %PBCH DM-RS MAPPING
        function addPbchDmRsToResourceGrid(obj,NCellId,pbchDmRs,t_offset,f_offset)
            arguments
                obj
                NCellId
                pbchDmRs
                t_offset = 0
                f_offset = 0
            end
            % nu parameter for shift of DM-RS
            nu=cast(mod(NCellId,4),"double");
            % pre-mapping staff
            pbchDmRs=obj.preparePbchDmRs(pbchDmRs);
            % first 60 PBCH DM-RS
            dmrs=pbchDmRs(1:60);
            
            % indexes array
            indexes=find(mod(1:1:240,4)==1);
            % mapping 1st part
            obj.resourceGrid(indexes+nu+f_offset,1+t_offset)=dmrs;
            % d---d---d---d … d---d---
            
            % last dm-rs block
            dmrs=pbchDmRs(end-59:end);
            % mapping 2nd part
            obj.resourceGrid(indexes+nu+f_offset,3+t_offset)=dmrs;
            % d---d---d---d … d---d---

            % next dm-rs block (24 elements are splitted into two blocks)
            dmrs=pbchDmRs(62:62+23);
            indexes=indexes(indexes<46 | indexes>192); % throwing SSS area
            % mapping 3rd part
            obj.resourceGrid(indexes+nu+f_offset,2+t_offset)=dmrs;
            % d---d---d--…-SSS-…-d---d--d
        end
    end
    
end