from aiogram.types import ReplyKeyboardRemove, \
    ReplyKeyboardMarkup, KeyboardButton, \
    InlineKeyboardMarkup, InlineKeyboardButton

sort_garbage = 'Сортировка_мусора'
habbit = "Привычки"

button_sort = KeyboardButton("/"+sort_garbage) #Текст приветствия, у кнопки обязательный параметр -- текст
button_habbits = KeyboardButton("/"+habbit)

main_kb = ReplyKeyboardMarkup(resize_keyboard=True, one_time_keyboard=True) #разметка клавы для ответа
main_kb.row(button_habbits, button_sort)
#main_kb.add(button_habbits)
#main_kb.add(button_sort)

