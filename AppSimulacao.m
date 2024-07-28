% Interface gráfica para simulação do controlador DMRNAC
% Autor: Christian Danner Ramos de Carvalho.

classdef AppSimulacao < matlab.apps.AppBase

    % Propriedades que correspondem aos componentes UI
    properties (Access = private)
        UIFigure            matlab.ui.Figure
        MatrixPanel         matlab.ui.container.Panel
        AEditFieldLabel     matlab.ui.control.Label
        AEditField          matlab.ui.control.EditField
        BEditFieldLabel     matlab.ui.control.Label
        BEditField          matlab.ui.control.EditField
        CEditFieldLabel     matlab.ui.control.Label
        CEditField          matlab.ui.control.EditField
        DEditFieldLabel     matlab.ui.control.Label
        DEditField          matlab.ui.control.EditField
        ParametersPanel     matlab.ui.container.Panel
        gammaEditFieldLabel matlab.ui.control.Label
        gammaEditField      matlab.ui.control.NumericEditField
        sigmaEditFieldLabel matlab.ui.control.Label
        sigmaEditField      matlab.ui.control.NumericEditField
        lambdaEditFieldLabel matlab.ui.control.Label
        lambdaEditField     matlab.ui.control.NumericEditField
        nEditFieldLabel     matlab.ui.control.Label
        nEditField          matlab.ui.control.NumericEditField
        bEditFieldLabel     matlab.ui.control.Label
        bEditField          matlab.ui.control.NumericEditField
        ftEditFieldLabel    matlab.ui.control.Label
        ftEditField         matlab.ui.control.NumericEditField
        dtEditFieldLabel    matlab.ui.control.Label
        dtEditField         matlab.ui.control.NumericEditField
        UncertaintyPanel    matlab.ui.container.Panel
        DeltaEditFieldLabel matlab.ui.control.Label
        DeltaEditField      matlab.ui.control.EditField
        Button              matlab.ui.control.Button
        Button2             matlab.ui.control.Button
        PlotAxes            matlab.ui.control.UIAxes
        PlotAxes2           matlab.ui.control.UIAxes
        PlotAxes3           matlab.ui.control.UIAxes
        PlotAxes4           matlab.ui.control.UIAxes
        LineSwitch          matlab.ui.control.Switch
    end

    % Métodos que compõem a UI
    methods (Access = private)

        % Código que é executado após a criação do objeto
        function startupFcn(app)
            % Configuração da janela principal
            app.UIFigure = uifigure('Name', 'SIMULADOR DMRNAC', 'Position', [100, 100, 820, 620]);

            % Painel para as matrizes do sistema dinâmico
            app.MatrixPanel = uipanel(app.UIFigure);
            app.MatrixPanel.Title = 'Dynamic System Matrices';
            app.MatrixPanel.Position = [50, 420, 350, 170]; % [margem esq, margem inferior, largura, altura]

            app.AEditFieldLabel = uilabel(app.MatrixPanel);
            app.AEditFieldLabel.Position = [10, 120, 100, 22];
            app.AEditFieldLabel.Text = 'A:';

            app.AEditField = uieditfield(app.MatrixPanel, 'text');
            app.AEditField.Position = [120, 120, 200, 22];
            app.AEditField.Value = '0, 1; 0, -0.1190';

            app.BEditFieldLabel = uilabel(app.MatrixPanel);
            app.BEditFieldLabel.Position = [10, 90, 100, 22];
            app.BEditFieldLabel.Text = 'B:';

            app.BEditField = uieditfield(app.MatrixPanel, 'text');
            app.BEditField.Position = [120, 90, 200, 22];
            app.BEditField.Value = '0; 373.3151';

            app.CEditFieldLabel = uilabel(app.MatrixPanel);
            app.CEditFieldLabel.Position = [10, 60, 100, 22];
            app.CEditFieldLabel.Text = 'C:';

            app.CEditField = uieditfield(app.MatrixPanel, 'text');
            app.CEditField.Position = [120, 60, 200, 22];
            app.CEditField.Value = '1, 0';

            app.DEditFieldLabel = uilabel(app.MatrixPanel);
            app.DEditFieldLabel.Position = [10, 30, 100, 22];
            app.DEditFieldLabel.Text = 'D:';

            app.DEditField = uieditfield(app.MatrixPanel, 'text');
            app.DEditField.Position = [120, 30, 200, 22];
            app.DEditField.Value = '0';         

            % Painel para os parâmetros editáveis
            
            app.ParametersPanel = uipanel(app.UIFigure);
            app.ParametersPanel.Title = 'Parameters';
            app.ParametersPanel.Position = [50, 125, 250, 284];

            app.gammaEditFieldLabel = uilabel(app.ParametersPanel);
            app.gammaEditFieldLabel.Position = [10, 234, 100, 22];
            app.gammaEditFieldLabel.Text = 'gamma:';

            app.gammaEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.gammaEditField.Position = [120, 234, 100, 22];
            app.gammaEditField.Value = 2;

            app.sigmaEditFieldLabel = uilabel(app.ParametersPanel);
            app.sigmaEditFieldLabel.Position = [10, 204, 100, 22];
            app.sigmaEditFieldLabel.Text = 'sigma:';

            app.sigmaEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.sigmaEditField.Position = [120, 204, 100, 22];
            app.sigmaEditField.Value = 500;

            app.lambdaEditFieldLabel = uilabel(app.ParametersPanel);
            app.lambdaEditFieldLabel.Position = [10, 174, 100, 22];
            app.lambdaEditFieldLabel.Text = 'lambda:';

            app.lambdaEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.lambdaEditField.Position = [120, 174, 100, 22];
            app.lambdaEditField.Value = 0.1;
          
            app.nEditFieldLabel = uilabel(app.ParametersPanel);
            app.nEditFieldLabel.Position = [10, 144, 100, 22];
            app.nEditFieldLabel.Text = 'n:';

            app.nEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.nEditField.Position = [120, 144, 100, 22];
            app.nEditField.Value = 500;

            app.bEditFieldLabel = uilabel(app.ParametersPanel);
            app.bEditFieldLabel.Position = [10, 114, 100, 22];
            app.bEditFieldLabel.Text = 'b:';

            app.bEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.bEditField.Position = [120, 114, 100, 22];
            app.bEditField.Value = 100;

            app.dtEditFieldLabel = uilabel(app.ParametersPanel);
            app.dtEditFieldLabel.Position = [10, 84, 100, 22];
            app.dtEditFieldLabel.Text = 'dt:';

            app.dtEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.dtEditField.Position = [120, 84, 100, 22];
            app.dtEditField.Value = 0.001;

            app.ftEditFieldLabel = uilabel(app.ParametersPanel);
            app.ftEditFieldLabel.Position = [10, 54, 100, 22];
            app.ftEditFieldLabel.Text = 'ft:';

            app.ftEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.ftEditField.Position = [120, 54, 100, 22];
            app.ftEditField.Value = 20;
          
            % Painel para a incerteza editável
            
            app.UncertaintyPanel = uipanel(app.UIFigure);
            app.UncertaintyPanel.Title = 'Uncertainty';
            app.UncertaintyPanel.Position = [50, 50, 350, 60];

            app.DeltaEditFieldLabel = uilabel(app.UncertaintyPanel);
            app.DeltaEditFieldLabel.Position = [10, 8, 100, 22];
            app.DeltaEditFieldLabel.Text = 'Delta:';

            app.DeltaEditField = uieditfield(app.UncertaintyPanel, 'text');
            app.DeltaEditField.Position = [50, 8, 290, 22];
            app.DeltaEditField.Value = '1+x(1)+x(2)+x(1)^2+sin(x(1))+cos(x(1))+sin(x(2))+cos(x(2))';

            % Botão para executar a simulação
            app.Button = uibutton(app.UIFigure, 'push');
            app.Button.Text = 'Run Simulation';
            app.Button.Position = [50, 10, 130, 30];
            app.Button.ButtonPushedFcn = @(btn,event) plotarGrafico(app);

            % Botão 2 para plotar os gráficos
            app.Button2 = uibutton(app.UIFigure, 'push');
            app.Button2.Text = 'Clear Chart';
            app.Button2.Position = [190, 10, 110, 30];
            app.Button2.ButtonPushedFcn = @(btn,event) limparGrafico(app);

            % Criar o switch
            app.LineSwitch = uiswitch(app.UIFigure, 'slider');
            app.LineSwitch.Position = [330, 15, 45, 30];
            app.LineSwitch.Items = {'Off', 'On'};
            app.LineSwitch.ValueChangedFcn = @(swi,event) LineSwitchValueChanged(app);

            % Área para plotar os gráficos
            app.PlotAxes = uiaxes(app.UIFigure);
            app.PlotAxes.Position = [400, 450, 400, 150];

            app.PlotAxes2 = uiaxes(app.UIFigure);
            app.PlotAxes2.Position = [400, 300, 400, 150];

            app.PlotAxes3 = uiaxes(app.UIFigure);
            app.PlotAxes3.Position = [400, 150, 400, 150];

            app.PlotAxes4 = uiaxes(app.UIFigure);
            app.PlotAxes4.Position = [400, 0, 400, 150];

        end

        % Função chamada quando o botão é pressionado
        function [t_rec,r_rec,x_rec,xmi_rec,u_rec,delta_rec,w_theta_rec] = runSimulation(app)
            % Obter os valores das entradas da matriz

            A = str2num(app.AEditField.Value);
            B = str2num(app.BEditField.Value);
            C = str2num(app.CEditField.Value);
            D = str2num(app.DEditField.Value);

            % Obter os valores dos parâmetros

            gamma = app.gammaEditField.Value;
            sigma = app.sigmaEditField.Value;
            lambda = app.lambdaEditField.Value;
            n = app.nEditField.Value;
            b = app.bEditField.Value;
            dt = app.dtEditField.Value;
            ft = app.ftEditField.Value;       

            % Condições iniciais
            x = [0; 0];
            xm = [0; 0];
            xmi = [0; 0];
            Psi = [0; 0]; % Vetor de estado do governador de comando
            
            % Controlador Nominal
            % K1 = lqr(A,B,eye(2),0.1);
            K1 = place(A,B,[-1 -2]); % Ganho de realimentação de estados
            K2 = -inv(C*inv(A-B*K1)*B); % Ganho de alimentação direta da referência

            % Condições de Correspondência do Modelo
            Am = A-B*K1;
            Bm = B*K2;

            % Sinal de correção (Recuperação de Desempenho)
            G = inv(K2)*inv(B'*B)*B'; % Matriz do sinal governador de comando
            Omega = B*inv(B'*B)*B';

            % Equação algébrica de Lyapunov
            Q = eye(2);
            P = lyap(Am',Q); %solução da equação algébrica de Lyapunov

            % Condições iniciais para W_hat
            W_hat = zeros(2*n+1,1);

            % Parâmetros da RNA
            centers = linspace(-b, b, n); % Vetor de centros dos neurônios
            Theta = zeros(2*n+1, 1); % Inicialização do vetor regressor da RNA

            % msgbox(['Texto digitado: ', A], 'Texto na Caixa');

            % Executa a iteração
            index = 1;
            for k = 0:dt:ft
                Delta = eval(app.DeltaEditField.Value);
                % Referência
                % if k <= 10
                %     r = 1;
                % end
                % if k > 10
                %     r = -1;
                % end
                if k < 2
                    r = 0;
                end
                if k >= 2
                    r = 1;
                end
                if k >= 6
                    r = -1;
                end
                if k >= 10
                    r = 1;
                end
                if k >= 14
                    r = -1;
                end

                % Montagem do vetor regressor da RNA
                for i = 1:n
                    Theta(i) = exp(-0.25*(abs(x(1)-centers(i)))^2); % RBFs
                    Theta(i+n) = exp(-0.25*(abs(x(2)-centers(i)))^2);
                end
                Theta(end) = 1; % Bias
                % Theta(1) = 1; % Bias
                % for i = 1:n
                %     Theta(i+1) = exp(-0.25*(abs(x(1)-centers(i)))^2); % RBFs
                %     Theta(i+n+1) = exp(-0.25*(abs(x(2)-centers(i)))^2);
                % end     

                % Sinal de correção para recuperação
                Psi = Psi + dt*(-lambda*(Psi-(x - xm))); % Dinâmica do governador de comando
                v = lambda*Psi + (Am-lambda*eye(2))*((x-xm)); % Sinal do governador de comando

                % Sinal de controle
                u = -K1*x + K2*(r+G*v) - W_hat'*Theta;
                xm = xm + dt*(Am*xm + Bm*r + Omega*v); % modelo de referência
                xmi = xmi + dt * (Am*xmi + Bm*r); % apenas para plotagem do modelo ideal
                W_hat = W_hat + dt*(gamma*(Theta*(x-xm)'*P*B - sigma*W_hat)); % Atualização dos pesos (adaptação)

                % Sistema Atual
                x = x + dt*(A*x + B*(u + Delta));

                % Gravação dos dados
                delta_rec(index,1) = Delta;
                w_theta_rec(index,1) = W_hat'*Theta;
                r_rec(index,1) = r;
                xm_rec(index,1:2) = xm;
                xmi_rec(index, 1:2) = xmi;
                x_rec(index,1:2) = x;
                u_rec(index,1) = u;
                t_rec(index,1) = k;
                e_rec(index,1:2) = x-xmi;
                w_til_rec(index,1) = Delta - W_hat'*Theta;
                index = index + 1;
            end

            return

        end

        % Função de callback para o switch
        function [lineColor] = LineSwitchValueChanged(app)

                % Obter o estado atual do switch
                switchState = app.LineSwitch.Value;

                % Determinar a cor da linha com base no estado do switch
                if switchState == "On"
                    lineColor = ''; % Outra cor
                else
                    lineColor = 'k-'; % Cor preta
                end

                return
        end

        % Função chamada quando o botão é pressionado
        function plotarGrafico(app)
            [lineColor] = LineSwitchValueChanged(app);
            [t_rec,r_rec,x_rec,xmi_rec,u_rec,delta_rec,w_theta_rec] = runSimulation(app);
            % Plota os gráficos atualizados
            hold(app.PlotAxes, 'on'); box(app.PlotAxes, 'on'); grid(app.PlotAxes, 'on');
            % title(app.PlotAxes,'DMNAC, mod-$\sigma$, sinal GC','fontsize',16,'interpreter','latex');
            p0 = plot(app.PlotAxes,t_rec,r_rec,'r:');
            set(p0,'linewidth',4);
            p1 = plot(app.PlotAxes,t_rec, xmi_rec(:,1), 'c-');
            set(p1, 'linewidth', 3);
            p2 = plot(app.PlotAxes,t_rec,x_rec(:,1),lineColor);
            set(p2,'linewidth',2);
            %text(app.PlotAxes, 5, max(x_rec(:,1)),'lambda = 10', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
            xlabel(app.PlotAxes,'$t$ (s)','fontsize',10,'interpreter','latex');
            ylabel(app.PlotAxes,'$x_1$ (m)','fontsize',10,'interpreter','latex');
            legend(app.PlotAxes,'$r$','$x(1)_m$','x(1)','fontsize',10,'interpreter','latex');
                        
            hold(app.PlotAxes2, 'on'); box(app.PlotAxes2, 'on'); grid(app.PlotAxes2, 'on');
            p1 = plot(app.PlotAxes2,t_rec,xmi_rec(:,2),'c-');
            set(p1,'linewidth',3);
            p2 = plot(app.PlotAxes2,t_rec,x_rec(:,2),lineColor);
            set(p2,'linewidth',2);
            xlabel(app.PlotAxes2,'$t$ (s)','fontsize',10,'interpreter','latex');
            ylabel(app.PlotAxes2,'$x_2$ (m/s)','fontsize',10,'interpreter','latex');
            legend(app.PlotAxes2,'$x(2)_m$','$x(2)$','fontsize',10,'interpreter','latex');

            hold(app.PlotAxes3, 'on'); box(app.PlotAxes3, 'on'); grid(app.PlotAxes3, 'on');
            p3 = plot(app.PlotAxes3,t_rec,u_rec,'b-');
            set(p3,'linewidth',3);
            xlabel(app.PlotAxes3,'$t$ (s)','fontsize',10,'interpreter','latex');
            ylabel(app.PlotAxes3,'$u(t)$','fontsize',10,'interpreter','latex');
            legend(app.PlotAxes3,'$u(t)$','fontsize',10,'interpreter','latex');

            hold(app.PlotAxes4, 'on'); box(app.PlotAxes4, 'on'); grid(app.PlotAxes4, 'on');
            plot(app.PlotAxes4,t_rec,delta_rec,'g-','linewidth',3); %hold on; box on; grid;
            plot(app.PlotAxes4,t_rec,w_theta_rec,'r--','linewidth',3); %grid;
            xlabel(app.PlotAxes4,'$t$ (s)','fontsize',10,'interpreter','latex');
            ylabel(app.PlotAxes4,'$\Delta(x)$','fontsize',10,'interpreter','latex');
            legend(app.PlotAxes4,'$\Delta(x)$','$\hat{W}^T\Theta(x)$','fontsize',10,'interpreter','latex');

            % Plotar os gráficos na área de plotagem
            % Dados para o gráfico
            % x = linspace(0, 10, 100);
            % y = sin(x);
            % 
            % % Plotar o gráfico
            % plot(app.PlotAxes, x, y);
            % grid(app.PlotAxes, 'on');
            % title(app.PlotAxes, 'Gráfico de y = sin(x)');
            % xlabel(app.PlotAxes, 'x');
            % ylabel(app.PlotAxes, 'y');
        end

        function limparGrafico(app)
            cla(app.PlotAxes); % Limpar o gráfico antes de traçar novos dados
            cla(app.PlotAxes2);
            cla(app.PlotAxes3);
            cla(app.PlotAxes4);
        end

    end

    % Métodos que são chamados quando certas ações ocorrem
    methods (Access = public)

        % Construção do app
        function app = AppSimulacao
            % Criação e configuração de componentes UI
            startupFcn(app);

            % Exibir a janela principal
            app.UIFigure.Visible = 'on';
        end
    end
end
