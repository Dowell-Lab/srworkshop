# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH

export LS_COLORS="$LS_COLORS:di=1;36"

export PS1="\[\e[0;33m\]\u@srAWS\[\e[0;37m\]:\[\e[1;36m\]\w \[\e[00m\]\$ "
