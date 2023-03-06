"""
@author: Mad Dreamer
"""

from tkinter import *
from PIL import ImageTk, Image

# Labels
LABEL_FONT_NAME = "Verdana"
LABEL_FONT_SIZE = 10
# Frames
FRAME_PADX = 1
FRAME_PADY = 1
FRAME_IPADX = 2
FRAME_IPADY = 1
# Buttons
BUTTON_FONT_NAME = "Verdana"
BUTTON_FONT_SIZE = 10
# Entry boxes
ENTRY_FONT_NAME = "Verdana"
ENTRY_FONT_SIZE = 10
# Check boxes
CHECK_FONT_NAME = "Verdana"
CHECK_FONT_SIZE = 10
# Check boxes
RADIO_FONT_NAME = "Verdana"
RADIO_FONT_SIZE = 10

root = Tk()
root.title("ME362")
root.geometry("1500x900")
#root.state("zoomed") # fullscreen
w = 1500 # width for the Tk root
h = 900 # height for the Tk root

# get screen width and height
ws = root.winfo_screenwidth() # width of the screen
hs = root.winfo_screenheight() # height of the screen

# calculate x and y coordinates for the Tk root window
#x = (ws/2) - (w/2)
#y = (hs/2) - (h/2)
x = 0
y = 0

# set the dimensions of the screen 
# and where it is placed
root.geometry('%dx%d+%d+%d' % (w, h, x, y))

###############################################################################
###############################################################################
###############################################################################
###############################################################################

test_label = Label(root, text="This is a label", font=(LABEL_FONT_NAME,LABEL_FONT_SIZE)).grid(row=0, column=0)
test_label2 = Label(root, text="This is a second label", font=(LABEL_FONT_NAME,LABEL_FONT_SIZE)).grid(row=2, column=0)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

test_frame = LabelFrame(root)
test_frame.grid(row=0, column=1, padx=FRAME_PADX, pady=FRAME_PADY, ipadx=FRAME_IPADX, ipady=FRAME_IPADY)
test_frame_label = Label(test_frame, text="This is a label in frame", font=(LABEL_FONT_NAME,LABEL_FONT_SIZE)).grid(row=0,column=0)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def test_button_action(n):
    print('this button works {}.'.format(n))

test_button = Button(root, text="My Button", font=(BUTTON_FONT_NAME, BUTTON_FONT_SIZE),
                     command=lambda:test_button_action(entry.get()),
                     #fg="black", bg="grey", activebackground="red"
                     )
test_button.grid(row=3,column=1)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

entry = Entry(root, width=20, font=(ENTRY_FONT_NAME,ENTRY_FONT_SIZE), justify="center")
entry.grid(row=0,column=2)

entry2 = Entry(root, width=20, font=(ENTRY_FONT_NAME,ENTRY_FONT_SIZE), justify="center",state="readonly")
entry2.grid(row=1,column=2)
entry2.configure(state="normal")
entry2.insert(0, "23")
entry2.configure(state="readonly")

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def check_action(var):
    print(var)

varCheck = StringVar()
varCheck.set("On")
checkbox = Checkbutton(root, font=(CHECK_FONT_NAME,CHECK_FONT_SIZE), text="Switch", onvalue="On", offvalue="Off", variable=varCheck,
                       command=lambda: check_action(varCheck.get())
                       )
checkbox.grid(row=3,column=0)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def scaler(val):
    global slider2
    global sliderentry
    sliderentry.configure(state="normal")
    sliderentry.delete(0, END)
    sliderentry.insert(0, f"{10**(slider2.get()/10):.2f}")
    sliderentry.configure(state="readonly")

sliderFrame = LabelFrame(root)
sliderFrame.grid(row=4, column=0, padx=FRAME_PADX, pady=FRAME_PADY, ipadx=FRAME_IPADX, ipady=FRAME_IPADY)

slider = Scale(sliderFrame, from_=-20, to=20, orient="horizontal",
               width=20, length=200)
slider.grid(row=0,column=0)

slider2 = Scale(sliderFrame, from_=-20, to=20, orient="horizontal",
               width=20, length=200,
               showvalue=0,
               command=scaler)
slider2.grid(row=1,column=0)

sliderentry = Entry(sliderFrame, width=5, font=(ENTRY_FONT_NAME,ENTRY_FONT_SIZE), justify="center",state="readonly")
sliderentry.grid(row=2,column=0)
sliderentry.configure(state="normal")
sliderentry.insert(0, 1.00)
sliderentry.configure(state="readonly")

###############################################################################
###############################################################################
###############################################################################
###############################################################################

radioFrame = LabelFrame(root)
radioFrame.grid(row=5, column=0, padx=FRAME_PADX, pady=FRAME_PADY, ipadx=FRAME_IPADX, ipady=FRAME_IPADY)

colorMapVar = StringVar()
colorMapVar.set("jet")

#imageJet = Image.open("jet.png")
#imageJet = ImageTk.PhotoImage(imageJet.resize((100,20)))
#imageCopper = Image.open("copper.png")
#imageCopper = ImageTk.PhotoImage(imageCopper.resize((100,20)))
#imageCool = Image.open("cool.png")
#imageCool = ImageTk.PhotoImage(imageCool.resize((100,20)))

radioJet = Radiobutton(radioFrame, variable=colorMapVar, text="Jet", value="jet",
                       font=(RADIO_FONT_NAME,RADIO_FONT_SIZE)
                       )
radioJet.grid(row=0, column=0, sticky="W")
#radioJet.configure(image=imageJet)

radioCopper = Radiobutton(radioFrame, variable=colorMapVar, text="Copper", value="copper",
                       font=(RADIO_FONT_NAME,RADIO_FONT_SIZE)
                       )
radioCopper.grid(row=1, column=0, sticky="W")
#radioCopper.configure(image=imageCopper)

radioCool = Radiobutton(radioFrame, variable=colorMapVar, text="Cool", value="cool",
                       font=(RADIO_FONT_NAME,RADIO_FONT_SIZE)
                       )
radioCool.grid(row=2, column=0, sticky="W")
#radioCool.configure(image=imageCool)


###############################################################################

root.mainloop()