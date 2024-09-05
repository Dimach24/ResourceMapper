classdef ResourceMapper<handle
    
    properties
        resourceGrid
        % main grid
        mu = 1
        % subcarrier spacing cofiguration (Δf=2^mu*15 [kHz]) that equals "1" for case 'C'
        frameCount
        % amount of frames in resourceGrid
    end
    
    methods
        
        function obj=createResourceGrid(obj, ...
                frameCount, ...
                channelBandwidth)
            % createResourceGrid
            % creates empty Resource grid for this
            % configuration [38.211, 4.3.2] or wipes
            % all data in the grid
            arguments
                obj                     ResourceMapper
                frameCount              (1,1)
                % amounts of empty frames to create
                channelBandwidth (1,1) = 60
                % channel bandwidth, MHz see [38.101-1: Table 5.3.2-1]
            end
            R_GRID_CONSTANTS;
            scs = 30; % subcarrier spacing in kHz
            
            NRB_tran_band_seq=MaximumTransmissionBandwidthConfiguration(scs);
            NRB=NRB_tran_band_seq(channelBandwidth);
            obj.resourceGrid=zeros(12*NRB,2^obj.mu*140*frameCount); % FIXME
            % initializing
            obj.frameCount=frameCount;
        end
        
        %%% ===============================================================
        function addSsBlockByCase(obj,nCellId,pssSignal,sssSignal,pbch,pbchDmRs,t_offset,f_offset,beta)
            % addSsBlockByCase
            % places signals according to  [38.213, 4.1]
            arguments
                obj ResourceMapper
                nCellId
                % physical layer cell identity; see [38.211,7.4.2.1]
                pssSignal (1,127)
                % array of pss values [38.211,7.4.2.2]
                sssSignal (1,127)
                % array of sss calues [38.211,7.4.2.2]
                pbch
                % physical broadcast channel data
                pbchDmRs
                % pbch demodulation reference signal [38.211,7.4.1.4]
                t_offset (1,1)
                % time domain offset
                f_offset (1,1)
                % freq. domain offset
                beta (1,4) = [1 1 1 1]
                % power allocation scaling factor
            end
            shifts = reshape([2,8]+14*[0 1 2 3].',1,[]); % SS/BCH block indexes for current case 'C'
            indexInData=1;
            % for each half-frame
            halfShifts=0:(2*obj.frameCount-1);
            halfShifts=halfShifts*2^obj.mu*70;
            for halfFrameShift=halfShifts
                for shift=(shifts+halfFrameShift)
                    addSsBlockToResourceGrid(obj,nCellId,pssSignal,sssSignal,pbch(indexInData,:),pbchDmRs(indexInData,:),t_offset+shift,f_offset,beta)
                    indexInData=indexInData+1;
                end
            end
        end
        
        function addSsBlockToResourceGrid(obj,nCellId,pssSignal,sssSignal,pbch,pbchDmRs,t_offset,f_offset,beta)
            % addSsBlockToResourceGrid
            % places signals according to configuration [38.211, 7.4.3.1-1]
            arguments
                obj ResourceMapper
                nCellId int32
                % physical layer cell identity; see [38.211,7.4.2.1]
                pssSignal
                % primary sync. signal [38.211, 7.4.2.2]
                sssSignal
                % secondary sync. signal [38.211, 7.4.2.3]
                pbch
                % physical broadcast channel data
                pbchDmRs
                % pbch demodulation reference signal [38.211,7.4.1.4]
                t_offset=0
                % time domain offset
                f_offset=0
                % freq. domain offset
                beta (1,4) = [1 1 1 1]
                % power allocation scaling factor
            end
            obj.addPssToResourceGrid(pssSignal,t_offset,f_offset,beta(1));
            obj.addSssToResourceGrid(sssSignal,t_offset+2,f_offset,beta(2));
            obj.addPbchToResourceGrid(nCellId,pbch,t_offset+1,f_offset,beta(3));
            obj.addPbchDmRsToResourceGrid(nCellId,pbchDmRs,t_offset+1,f_offset,beta(4));
        end
        
        % SS MAPPING
        function addPssToResourceGrid(obj,PssSignal,t_offset,f_offset,beta)
            % adds PSS to resource matrix
            arguments
                obj
                PssSignal (1,127)
                % primary sync. signal [38.211, 7.4.2.2]
                t_offset (1,1) = 0
                % time domain offset
                f_offset (1,1) = 0
                % freq. domain offset
                beta (1,1) = 1
                % power allocation factor
            end
            obj.resourceGrid((57:183)+f_offset,1+t_offset) = beta .* PssSignal;
        end
        
        function  addSssToResourceGrid(obj, SssSignal,t_offset,f_offset,beta)
            % adds SSS to resource matrix
            arguments
                obj
                SssSignal (1,127)
                % secondary sync. signal [38.211, 7.4.2.3]
                t_offset (1,1) = 0
                % time domain offset
                f_offset (1,1) = 0
                % freq. domain offset
                beta (1,1) = 1
                % power allocation factor
            end
            obj.resourceGrid((57:183)+f_offset,1+t_offset) = beta .* SssSignal;
        end
        
        function addPbchToResourceGrid(obj,NCellId,pbch,t_offset,f_offset,beta)
            arguments
                obj ResourceMapper
                NCellId (1,1)
                % physical layer cell identity; see [38.211,7.4.2.1]
                pbch
                % physical broadcast channel data
                t_offset = 0
                % time domain offset
                f_offset = 0
                % freq. domain offset
                beta (1,1) = 1
                % power allocation factor
            end
            
            % nu parameter for shift of DM-RS
            nu=mod(NCellId,4);
            
            %throwing out dmrs indexes
            indexes=find(mod(0:1:239,4)~=(nu));
            
            % mapping first 180 PBCH
            obj.resourceGrid(indexes+f_offset,1+t_offset)=beta .* pbch(1:180);
            % mapping last 180
            obj.resourceGrid(indexes+f_offset,3+t_offset)=beta .* pbch(end-179:end);
            % mapping arround SSS
            indexes=indexes(indexes<49 | indexes>192);
            obj.resourceGrid(indexes+f_offset,2+t_offset)=beta .* pbch(181:181+71);
        end
        
        % PBCH DM-RS MAPPING
        function addPbchDmRsToResourceGrid(obj,NCellId,pbchDmRs,t_offset,f_offset,beta)
            arguments
                obj ResourceMapper
                NCellId (1,1)
                % physical layer cell identity; see [38.211,7.4.2.1]
                pbchDmRs
                % pbch demodulation reference signal [38.211,7.4.1.4]
                t_offset = 0
                % time domain offset
                f_offset = 0
                % freq. domain offset
                beta(1,1) = 1
                % power allocation factor
            end
            % nu parameter for shift of DM-RS
            nu=cast(mod(NCellId,4),"double");

            % first 60 PBCH DM-RS
            dmrs=pbchDmRs(1:60);
            
            % indexes array
            indexes=0:4:236;
            % mapping 1st part
            obj.resourceGrid(indexes+nu+f_offset+1,1+t_offset)=beta .* dmrs;
            % d---d---d---d … d---d---
            
            % last dm-rs block
            dmrs=pbchDmRs(end-59:end);
            % mapping 2nd part
            obj.resourceGrid(indexes+nu+f_offset+1,3+t_offset)=beta .* dmrs;
            % d---d---d---d … d---d---
            
            % next dm-rs block (24 elements are splitted into two blocks)
            dmrs=pbchDmRs(62:62+23);
            indexes=[0:4:44,192:4:236];
            % mapping 3rd part
            obj.resourceGrid(indexes+nu+f_offset+1,2+t_offset)=beta .* dmrs;
            % d---d---d--…-SSS-…-d---d--d
        end
    end
end