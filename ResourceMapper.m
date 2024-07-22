classdef ResourceMapper<handle
    
    properties
        resourceGrid
        % main grid
        mu
        % subcarrier spacing cofiguration (Δf=2^mu*15 [kHz])
        frameCount
        % amount of frames in resourceGrid
        isCycledPrefixExtended
        % cycled prefix flag
    end
    
    methods
        
        function obj=createResourceGrid(obj, ...
                fcase, ...
                pointA,...
                frameCount, ...
                isCycledPrefixExtended, ...
                tran_bandwidth)
            % createResourceGrid
            % creates empty Resource grid for this
            % configuration [38.211, 4.3.2] or wipes
            % all data in the grid
            arguments
                obj                     ResourceMapper
                fcase                      char
                % freq. config case letter, see [38.213, 4.1]
                pointA                  (1,1)
                % resource grid carrier frequency in GHz
                frameCount              (1,1)
                % amounts of empty frames to create
                isCycledPrefixExtended  (1,1) = false
                % extended cycled prefix
                tran_bandwidth (1,1)= 60
                % transmission bandwidth, MHz see [38.101-1: Table 5.3.2-1]
            end
            R_GRID_CONSTANTS;
            [mu,scs] = obj.caseDecipher(fcase,pointA);
            
            NRB_tran_band_seq=MaximumTransmissionBandwidthConfiguration(scs);
            NRB=NRB_tran_band_seq(tran_bandwidth);
            if(isCycledPrefixExtended && (mu==2)) % ONLY FOR MU==2
                %   12 sym per slot, 10*4 slot per frame = 480 symb per frame
                obj.resourceGrid=zeros(12*NRB,480*frameCount); % FIXME
            else
                obj.resourceGrid=zeros(12*NRB,2^mu*140*frameCount);% FIXME
            end
            % initializing
            obj.mu=mu;
            obj.isCycledPrefixExtended=isCycledPrefixExtended;
            obj.frameCount=frameCount;
        end
        
        %%% ===============================================================
        function addSsBlockByCase(obj,fcase,pointA,nCellId,pssSignal,sssSignal,pbch,pbchDmRs,t_offset,f_offset,beta)
            % addSsBlockByCase
            % places signals according to  [38.213, 4.1]
            arguments
                obj ResourceMapper
                fcase char
                % freq. config case, see [38.213, 4.1]
                pointA
                % resource grid carrier frequency in GHz
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
            [mu,~,shifts] = obj.caseDecipher(fcase,pointA);
            indexInData=1;
            % for each half-frame
            halfShifts=0:(2*obj.frameCount-1);
            if obj.isCycledPrefixExtended
                halfShifts=halfShifts*240;
            else
                halfShifts=halfShifts*2^obj.mu*70;
            end
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
            indexes=find(mod(1:1:240,4)==1);
            % mapping 1st part
            obj.resourceGrid(indexes+nu+f_offset,1+t_offset)=beta .* dmrs;
            % d---d---d---d … d---d---
            
            % last dm-rs block
            dmrs=pbchDmRs(end-59:end);
            % mapping 2nd part
            obj.resourceGrid(indexes+nu+f_offset,3+t_offset)=beta .* dmrs;
            % d---d---d---d … d---d---
            
            % next dm-rs block (24 elements are splitted into two blocks)
            dmrs=pbchDmRs(62:62+23);
            indexes=[(0:11)*4,(48:59)*4]+nu+1;
            % mapping 3rd part
            obj.resourceGrid(indexes+nu+f_offset,2+t_offset)=beta .* dmrs;
            % d---d---d--…-SSS-…-d---d--d
        end
        
        function [mu, scs, shifts, Lmax_] = caseDecipher(obj, fcase, f_carrier, isSpectrumAccessShared, isSpectrumOperationPaired) 
        % defines resource grid configuration via case identificator letter  
            arguments
                obj ResourceMapper
                fcase char % freq. config case letter, see [38.213, 4.1]
                f_carrier (1,1) % carrier frequency in GHz
                isSpectrumAccessShared logical = 0; % defines shared spectrum channel acces cases
                isSpectrumOperationPaired logical = 0; % defines paired spectrum operation cases
            end
            switch fcase
                case 'A'
                    if isSpectrumAccessShared
                        n = [0 1 2 3 4];
                        Lmax_ = 10;
                    else
                        n = [0 1];
                        Lmax_ = 4;
                        if (f_carrier > 3)
                            n = [0 1 2 3];
                            Lmax_ = 8;
                        end
                    end
                    shifts=reshape([2 8]+14*n.',1,[]);
                    mu = 0;
                case 'B'
                    n = 0;
                    Lmax_ = 4;
                    if f_carrier > 3
                        n = [0 1];
                        Lmax_ = 8;
                    end
                    shifts=reshape([4 8 16 20]+28*n.',1,[]);
                    mu = 1;
                case 'C'
                    if isSpectrumAccessShared
                        n = [0 1 2 3 4 5 6 7 8 9];
                        Lmax_ = 20;
                    else
                        n = [0 1];
                        Lmax_ = 4;
                        if ((f_carrier > 3)*isSpectrumOperationPaired||(f_carrier >= 1.88)*~isSpectrumOperationPaired)
                            n = [0 1 2 3];
                            Lmax_ = 8;
                        end
                    end
                    shifts=reshape([2,8]+14*n.',1,[]);
                    mu = 1;
                case 'D'
                    n = [0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18];
                    shifts=reshape([8,12,16,20]+28*n.',1,[]);
                    mu = 3;
                    Lmax_ = 64;
                case 'E'
                    n = [0 1 2 3 5 6 7 8];
                    shifts=reshape([(8:4:20),(32:4:44)]+56*n.',1,[]);
                    mu = 4;
                    Lmax_ = 64;
                case 'F'
                    n = 0:31;
                    shifts=reshape([2,9]+14*n.',1,[]);
                    mu = 5;
                    Lmax_ = 64;
                case 'G'
                    n = 0:31;
                    shifts=reshape([2,9]+14*n.',1,[]);
                    mu = 6;
                    Lmax_ = 64;
                otherwise
                    throw(MException("Freq. Case Error","Case must be uppercase letter A...G"));
            end
            shifts=sort(shifts); % indexing from 1
            scs = 2^mu * 15;
        end
    end
end