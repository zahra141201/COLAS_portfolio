'use strict';
/** @type {import('sequelize-cli').Migration} */
module.exports = {
  async up(queryInterface, Sequelize) {
    await queryInterface.createTable('Requests', {
      id: {
        allowNull: false,
        autoIncrement: true,
        primaryKey: true,
        type: Sequelize.INTEGER
      },
      pick_up_date: {
        type: Sequelize.DATE
      },
      comment: {
        type: Sequelize.STRING
      },
      state: {
        type: Sequelize.STRING
      },
      made_by: {
        type: Sequelize.STRING,
        references: { model: 'Users', key: 'email' },
        onUpdate: 'CASCADE',
        onDelete: 'CASCADE'
      },
      id_ingrediente: {
        type: Sequelize.INTEGER,
        references: { model: 'Ingredients', key: 'id' },
        onUpdate: 'CASCADE',
        onDelete: 'CASCADE'
      },
      createdAt: {
        allowNull: false,
        type: Sequelize.DATE
      },
      updatedAt: {
        allowNull: false,
        type: Sequelize.DATE
      }
    });
  },
  async down(queryInterface, Sequelize) {
    await queryInterface.dropTable('Requests');
  }
};
