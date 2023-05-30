classdef MixParams < handle

    % Sine-Mix parameters

    properties
        Fs
        ws
        dc
        a
        p
        ph

        has_explicit_dc
    end

    methods
        function obj = MixParams(params)

            arguments
                params.Fs
                params.ws
                params.dc = []
                params.a
                params.p = []
                params.ph = []
                params.v = []
            end

            Fs = params.Fs;
            ws = params.ws;
            a = params.a;
            p = params.p;
            dc = params.dc;
            ph = params.ph;
            
            assert(isempty(params.v) ~= isempty(params.dc), "Must specify exactly one of dc and v");

            obj.has_explicit_dc = ~isempty(dc);

            if iscell(a)
                a = cell2mat(a(:));
            end
            
            nSamples = size(a,2);

            if isempty(p)
                assert(~isempty(ph))
                p = repmat(ph(:),[1,nSamples]);
            end
            if iscell(p)
                p = cell2mat(p(:));
            end
            if ~obj.has_explicit_dc
                dc = zeros([1,nSamples]);
            end

            assert(all(size(a)==[length(ws),nSamples]));
            assert(all(size(p)==[length(ws),nSamples]));
            assert(numel(dc)==nSamples);
            assert(isempty(ph) || length(ph)==length(ws));

            obj.Fs = Fs;
            obj.ws = ws(:)';
            obj.dc = dc(:)';
            obj.a = a;
            obj.p = p;
            obj.ph = ph;

            if ~obj.has_explicit_dc
                obj.dc = params.v(:)' - obj.v;
            end

        end
        
        function smooth(obj, win_sz)
            obj.a = movavg(obj.a, win_sz);
            obj.p = movavg(obj.p, win_sz);
        end

        function tt = t(obj)
            tt = (0:(length(obj.dc)-1))/obj.Fs;
        end

        function c = comps(obj)
            t = obj.t;
            c = nan([length(obj.ws)+1, length(t)]);
            c(1,:) = obj.dc;
            for i = 1 : length(obj.ws)
                c(i+1, :) = obj.a(i,:).*sin(obj.ws(i)*t + obj.p(i,:));
            end
            assert(~any(isnan(c(:))));
        end

        function vv = v(obj)
            vv = sum(obj.comps,1);
        end

        function aa = amps(obj)
            aa = [obj.dc; obj.a];
        end

        function set_dc_from_v(obj, v)
            obj.dc = zeros([1,length(v)]);
            obj.has_explicit_dc = false;
            obj.dc = v - obj.v;
        end

    end
        

end